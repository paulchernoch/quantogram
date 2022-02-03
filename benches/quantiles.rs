#[macro_use]
extern crate bencher;
extern crate quantogram;

use bencher::Bencher;
use quantogram::Quantogram;
use quantiles::ckms::CKMS;
use zw_fast_quantile::UnboundEpsilonSummary;

// Helpers

/// Exercises Quantogram::median
fn many_medians(data: &Vec<f64>, perform_median: bool) {
    let mut q = Quantogram::new();
    for sample in data.iter() {
        q.add(*sample);
        if perform_median{
            let _median = q.median().unwrap();
        }
    }
}

/// Exercises CKMS::query(0.5)
fn many_medians_ckms(data: &Vec<f64>, perform_median: bool) {
    let mut ckms = CKMS::<f64>::new(0.01);
    for sample in data.iter() {
        ckms.insert(*sample);
        if perform_median{
            let _median = ckms.query(0.5).unwrap();
        }
    }
}

/// Exercises UnboundEpsilonSummary::query(0.5)
fn many_medians_zw(data: &Vec<f64>, perform_median: bool) {
    let mut zw = UnboundEpsilonSummary::new(0.01);
    for sample in data.iter() {
        // Need an ordered type - f64 is only partially ordered. 
        // Change floats to integers for test.
        zw.update((*sample * 1000.0).round() as i64);
        if perform_median{
            let _median = zw.query(0.5);
        }
    }
}

fn randomish(bottom: f64, top:f64, count: usize) -> Vec<f64> {
    let mut data : Vec<f64> = Vec::new();
    let range = top - bottom;
    let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0;
    for i in 0..count {
        // Pseudo random Weyl sequence based on Golden ratio.
        let pseudo_rand = ((i as f64) * phi) % 1.0_f64;
        // By squaring pseudo_rand we get a non-uniform distribution,
        // which is more interesting.
        let sample = bottom + pseudo_rand * pseudo_rand * range;
        data.push(sample);
    }
    data
}

// Benchmarks

// Quantogram

fn bench_median_200000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 200000);
    b.iter(|| {  many_medians(&data, true); });  
}
fn bench_median_100000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians(&data, true); });  
}
fn bench_median_050000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians(&data, true); });  
}
fn bench_median_010000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians(&data, true); });  
}
fn bench_median_001000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 1000);
    b.iter(|| {  many_medians(&data, true); });  
}

// CKMS

fn bench_median_200000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 200000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}
fn bench_median_100000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}
fn bench_median_050000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}
fn bench_median_010000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}
fn bench_median_001000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 1000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}

// ZW

fn bench_median_001000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 1000);
    b.iter(|| {  many_medians_zw(&data, true); });  
}
fn bench_median_010000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians_zw(&data, true); });  
}
fn bench_median_050000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians_zw(&data, true); });  
}

fn bench_median_100000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians_ckms(&data, true); });  
}

// Inserts

fn bench_insert_010000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians_ckms(&data, false); });  
}
fn bench_insert_050000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians_ckms(&data, false); });  
}
fn bench_insert_100000_ckms(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians_ckms(&data, false); });  
}

fn bench_insert_010000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians(&data, false); });  
}
fn bench_insert_050000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians(&data, false); });  
}
fn bench_insert_100000_qg(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians(&data, false); });  
}

fn bench_insert_010000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 10000);
    b.iter(|| {  many_medians_zw(&data, false); });  
}
fn bench_insert_050000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 50000);
    b.iter(|| {  many_medians_zw(&data, false); });  
}
fn bench_insert_100000_zw(b: &mut Bencher) {
    let data = randomish(-10000.0, 10000.0, 100000);
    b.iter(|| {  many_medians_zw(&data, false); });  
}


benchmark_group!(benches, 
    bench_median_001000_qg,   bench_median_010000_qg,   bench_median_050000_qg,  bench_median_100000_qg, bench_median_200000_qg,
    bench_median_001000_ckms, bench_median_010000_ckms, bench_median_050000_ckms,
    bench_median_001000_zw,   bench_median_010000_zw,   bench_median_050000_zw, bench_median_100000_zw,

    bench_insert_010000_qg,   bench_insert_010000_ckms, bench_insert_010000_zw,
    bench_insert_050000_qg,   bench_insert_050000_ckms, bench_insert_050000_zw,
    bench_insert_100000_qg,   bench_insert_100000_ckms, bench_insert_100000_zw
);
benchmark_main!(benches);
