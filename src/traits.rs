use std::fmt::{Debug, Display};
use std::ops::{Add, Div, Mul, Neg, Sub};
use rand::Rng;

pub trait Group<S>:

Clone+ Debug+ Display
+ Add<Output = Self>
+ Mul<S, Output = Self>  // Scalar multiplication: S * T -> T
+ PartialEq
{
    fn zero(d: usize) -> Self;
    fn random<R: Rng>(rng: &mut R, nonce: u16, d:usize) -> Self;
    fn sample_from<R: Rng>(rng: &mut R, sigma: f64, d:usize) -> Self;
    fn norm(&self) -> i32;
    fn dimension(&self) -> usize;
}

// Trait for the scalar field (e.g., integers mod p).
// Supports field operations: add, sub, mul, div (via invert), neg.
pub trait ScalarField:
Clone+ Copy+ Debug+ Display
+ Add<Output = Self>
+ Sub<Output = Self>
+ Mul<Output = Self>
+ Div<Output = Self>
+ Neg<Output = Self>
+ PartialEq
+ From<u32>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn invert(&self) -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
    fn norm(&self) -> i32;
}