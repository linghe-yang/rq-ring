use std::cmp::PartialEq;
use crystals_dilithium::poly::*;
use crystals_dilithium::params;
use rand::Rng;
use serde::de::Visitor;
use serde::ser::SerializeSeq;
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::ops::{Add, Mul, Sub};
use crystals_dilithium::poly::{Poly, reduce, caddq};
use num_traits::Zero;
use crate::traits::{Group, ScalarField};
use crate::zq::Zq;

// Rq represents an element in Z_Q[x] / (x^N + 1)
pub const Q: i32 = params::Q;
pub const N: i32 = params::N;
#[derive(Clone, Copy)]
pub struct Rq {
    pub poly: Poly,
}

impl Rq {
    pub fn sample_discrete_gaussian<R: Rng>(rng: &mut R, sigma: f64) -> Self {
        assert!(sigma > 0.0, "Sigma must be positive");
        // assert!(
        //     sigma <= 6.0,
        //     "Sigma too large"
        // );
        let tail = (12.0 * sigma).ceil() as i32;
        let mut poly = Poly { coeffs: [0i32; params::N as usize] };
        for coeff in poly.coeffs.iter_mut() {
            // rejection sampling
            loop {
                let x = rng.gen_range(-tail..=tail);
                let prob = ((-(x as f64).powi(2)) / (2.0 * sigma.powi(2))).exp();
                let u = rng.r#gen::<f64>(); //[0, 1)
                if u < prob {
                    *coeff = x;
                    break;
                }
            }
        }
        reduce(&mut poly);
        Rq { poly }
    }

    // Serialize the Rq (Poly coefficients) to a byte vector
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity((params::N * 4) as usize); // Each i32 is 4 bytes
        let mut poly = self.poly;
        reduce(&mut poly);
        caddq(&mut poly);
        for coeff in poly.coeffs.iter() {
            bytes.extend(coeff.to_le_bytes());
        }
        bytes
    }

    // Deserialize a byte slice into an Rq, if valid
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != (params::N * 4) as usize { // Each i32 is 4 bytes
            return None;
        }
        let mut coeffs = [0i32; params::N as usize];
        for i in 0..params::N as usize {
            let start = i * 4;
            if start + 4 <= bytes.len() {
                let coeff_bytes: &[u8] = &bytes[start..start + 4];
                coeffs[i] = i32::from_le_bytes(coeff_bytes.try_into().ok()?);
            } else {
                return None;
            }
        }
        Some(Rq {
            poly: Poly { coeffs }
        })
    }
}


// Manual implementation of PartialEq for Rq
impl PartialEq for Rq {
    fn eq(&self, other: &Rq) -> bool {
        let mut t1 = self.poly;
        let mut t2 = other.poly;
        reduce(&mut t1);
        caddq(&mut t1);
        reduce(&mut t2);
        caddq(&mut t2);
        t1.coeffs == t2.coeffs
    }
}

impl PartialEq<Rq> for &Rq {
    fn eq(&self, other: &Rq) -> bool {
        let mut t1 = self.poly;
        let mut t2 = other.poly;
        reduce(&mut t1);
        caddq(&mut t1);
        reduce(&mut t2);
        caddq(&mut t2);
        t1.coeffs == t2.coeffs
    }
}

impl Zero for Rq {
    fn zero() -> Self {
        Rq {
            poly: Poly::default(),
        }
    }

    fn is_zero(&self) -> bool {
        self == <Rq as Group<Zq>>::zero(0)
    }
}

// Implementation of Group trait for Rq with Zq as scalar
impl Group<Zq> for Rq {
    // Zero element: polynomial with all coefficients 0
    fn zero(_: usize) -> Self {
        Rq {
            poly: Poly::default(),
        }
    }

    // Random element: generate polynomial with coefficients in [0, Q-1] using poly::uniform with provided nonce
    fn random<R: Rng>(rng: &mut R, nonce: u16, _:usize) -> Self {
        let mut poly = Poly::default();
        // Generate a random seed for uniform sampling
        let mut seed = [0u8; 32];
        rng.fill_bytes(&mut seed);
        uniform(&mut poly, &seed, nonce);
        Rq {
            poly,
        }
    }
    fn sample_from<R: Rng>(rng: &mut R, sigma: f64, _: usize) -> Self {
        Self::sample_discrete_gaussian(rng, sigma)
    }


    fn norm(&self) -> i32 {
        let mut max_norm = 0;
        let poly = self.poly.clone();

        // reduce(&mut poly);
        for &coeff in poly.coeffs.iter() {
            let mut t = coeff >> 31;
            t = coeff - (t & 2 * coeff);
            // Update maximum norm
            if t > max_norm {
                max_norm = t;
            }
        }
        max_norm
    }

    fn dimension(&self) -> usize {
        params::N as usize
    }
}

// Addition in Zq[x] / (x^n + 1): coefficient-wise addition with lazy reduction
impl Add for Rq {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let result = add(&self.poly, &other.poly);
        Rq {
            poly: result,
        }
    }
}

impl Sub for Rq {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        let result = sub(&self.poly, &other.poly);
        Rq {
            poly: result,
        }
    }
}

// Scalar multiplication: multiply each coefficient by a scalar in Zq
impl Mul<Zq> for Rq {
    type Output = Self;
    fn mul(self, scalar: Zq) -> Self {
        let mut result = Poly::default();
        let s: i64 = scalar.norm() as i64;
        for (i,x) in self.poly.coeffs.iter().enumerate() {
            let prod = *x as i64 * s;
            result.coeffs[i] = prod as i32;
        }
        reduce(&mut result);
        // for i in 0..params::N as usize {
        //     let prod = (self.poly.coeffs[i] as i64 * s) % params::Q as i64;
        //     result.coeffs[i] = prod as i32;
        // }
        Rq {
            poly: result,
        }
    }
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    let mut a_ntt = *a;
    let mut b_ntt = *b;
    ntt(&mut a_ntt);
    ntt(&mut b_ntt);
    let mut c_ntt = Poly::default();
    pointwise_montgomery(&mut c_ntt, &a_ntt, &b_ntt);
    invntt_tomont(&mut c_ntt);
    reduce(&mut c_ntt);
    c_ntt
}

// Multiplication in Zq[x] / (x^n + 1): polynomial multiplication using NTT
impl Mul for Rq {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let result = poly_mul(&self.poly, &other.poly);
        Rq {
            poly: result,
        }
    }
}
impl Mul<Rq> for Zq {
    type Output = Rq;
    fn mul(self, rq: Rq) -> Rq {
        rq * self // Reuse the existing impl
    }
}

// Manual implementation of Debug for Rq
impl fmt::Debug for Rq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rq")
            .field("coeffs", &self.poly.coeffs)
            .finish()
    }
}

// Implementation of Display for Rq (as polynomial)
impl fmt::Display for Rq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms = Vec::new();
        for (i, &coeff) in self.poly.coeffs.iter().enumerate() {
            if coeff != 0 {
                let coeff = Zq::new(coeff).norm();
                let sign = if coeff < 0 { "-" } else if i > 0 { "+" } else { "" };
                let c = coeff.to_string();
                let term = match i {
                    0 => format!("{}{}", sign, c),
                    1 => format!("{}{}x", sign, if coeff == 1 { "" } else { &*c }),
                    _ => format!("{}{}x^{}", sign, if coeff == 1 { "" } else { &*c }, i),
                };
                terms.push(term);
            }
        }
        if terms.is_empty() {
            write!(f, "0")
        } else {
            write!(f, "{}", terms.join(" "))
        }
    }
}


impl Serialize for Rq {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Use to_bytes to get the polynomial coefficients as a byte vector
        let bytes = self.to_bytes();
        // Serialize as a sequence of params::N * 4 bytes
        let mut seq = serializer.serialize_seq(Some(bytes.len()))?;
        for &byte in &bytes {
            seq.serialize_element(&byte)?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for Rq {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct RqVisitor;

        impl<'de> Visitor<'de> for RqVisitor {
            type Value = Rq;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str(&format!("a sequence of {} bytes representing an Rq polynomial", params::N * 4))
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: de::SeqAccess<'de>,
            {
                // Collect exactly params::N * 4 bytes
                let mut bytes = vec![0u8; (params::N * 4) as usize];
                for i in 0..bytes.len() {
                    bytes[i] = seq.next_element()?.ok_or_else(|| de::Error::invalid_length(i, &self))?;
                }
                // Use from_bytes to reconstruct Rq
                Rq::from_bytes(&bytes).ok_or_else(|| de::Error::custom("Invalid byte sequence: must represent a valid Rq polynomial"))
            }
        }

        deserializer.deserialize_seq(RqVisitor)
    }
}