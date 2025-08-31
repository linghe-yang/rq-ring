use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Mul, Sub, Div, Neg};
use rand::Rng;
use num_modular::{ModularInteger, MontgomeryInt};
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use serde::de::Visitor;
use serde::ser::SerializeTuple;
use crate::traits::ScalarField;

use super::rq::Q as DILITHIUM_Q;
const Q: u32 = DILITHIUM_Q as u32;
#[derive(Clone)]
pub struct Zq {
    pub(crate) value: MontgomeryInt<u32>,
}
impl Zq {
    pub fn new(val: i32) -> Self {
        let q_i32 = Q as i32;
        let mut v = val % q_i32;
        if v < 0 {
            v += q_i32;
        }
        Zq {
            value: MontgomeryInt::new(v as u32, &Q),
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.value.residue().to_le_bytes().to_vec()
    }

    // Deserialize a byte slice into a Zq, assuming it represents a reduced residue
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 4 { // u32 is 4 bytes
            return None;
        }
        let residue = u32::from_le_bytes(bytes.try_into().ok()?);
        if residue >= Q { // Ensure the residue is in [0, Q-1]
            return None;
        }
        Some(Zq {
            value: MontgomeryInt::new(residue, &Q),
        })
    }

}

impl ScalarField for Zq {
    fn zero() -> Self {
        Zq {
            value: MontgomeryInt::new(0, &Q),
        }
    }

    fn one() -> Self {
        Zq {
            value: MontgomeryInt::new(1, &Q),
        }
    }

    fn invert(&self) -> Self {
        Zq {
            value: self.value.inv().expect("Invert failed: input must be non-zero and coprime with Q"),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        let random_val = rng.gen_range(0..Q);
        Zq {
            value: MontgomeryInt::new(random_val, &Q),
        }
    }

    fn norm(&self) -> i32 {
        self.value.residue() as i32
    }

}

// 实现运算符 trait
impl Add for Zq {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Zq {
            value: self.value + other.value,
        }
    }
}

impl Sub for Zq {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Zq {
            value: self.value - other.value,
        }
    }
}

impl Mul for Zq {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Zq {
            value: self.value * other.value,
        }
    }
}

impl Div for Zq {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        Zq {
            value: self.value / other.value,
        }
    }
}

impl Neg for Zq {
    type Output = Self;

    fn neg(self) -> Self {
        Zq {
            value: self.value.neg(),
        }
    }
}

impl From<u32> for Zq {
    fn from(value: u32) -> Self {
        Zq {
            value: MontgomeryInt::new(value, &Q),
        }
    }
}

impl Copy for Zq {}

impl PartialEq for Zq {
    fn eq(&self, other: &Self) -> bool {
        self.value.residue() == other.value.residue()
    }
}

impl PartialOrd for Zq {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.value.residue().cmp(&other.value.residue()))
    }
}

impl Eq for Zq {}

impl fmt::Debug for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Zq({})", self.value.residue())
    }
}

impl fmt::Display for Zq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value.residue())
    }
}
impl Serialize for Zq {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // Use to_bytes to get the residue as a 4-byte vector
        let bytes = self.to_bytes();
        // Serialize as a tuple of 4 u8 values
        let mut tup = serializer.serialize_tuple(4)?;
        for &byte in &bytes {
            tup.serialize_element(&byte)?;
        }
        tup.end()
    }
}

impl<'de> Deserialize<'de> for Zq {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ZqVisitor;

        impl<'de> Visitor<'de> for ZqVisitor {
            type Value = Zq;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a 4-byte sequence representing a Zq residue")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: de::SeqAccess<'de>,
            {
                // Collect 4 bytes
                let mut bytes = [0u8; 4];
                for i in 0..4 {
                    bytes[i] = seq.next_element()?.ok_or_else(|| de::Error::invalid_length(i, &self))?;
                }
                // Use from_bytes to reconstruct Zq
                Zq::from_bytes(&bytes).ok_or_else(|| de::Error::custom("Invalid residue: must be in [0, Q-1]"))
            }
        }

        deserializer.deserialize_tuple(4, ZqVisitor)
    }
}