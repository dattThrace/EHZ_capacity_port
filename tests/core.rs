use ehz_capacity::*;
use nalgebra::{DMatrix, DVector};
#[test]
fn products_and_centroid() {
    assert_eq!(
        ktimes_t(&vec![vec![1.], vec![2.]], &vec![vec![3.], vec![4.]]),
        vec![vec![1., 3.], vec![2., 3.], vec![1., 4.], vec![2., 4.]]
    );
    let k = vec![vec![0., 0.], vec![2., 0.], vec![0., 3.]];
    let c = barycenter(&k, &[vec![0, 1, 2]]).unwrap();
    assert!((c[0] - 2. / 3.).abs() < 1e-12 && (c[1] - 1.).abs() < 1e-12);
    assert_eq!(omega(&[1., 0.], &[0., 1.]), 1.);
    assert_eq!(omega(&[1., 2., 3., 4.], &[5., 6., 7., 8.]), -16.);
}
#[test]
fn equality_branches() {
    let a = DMatrix::from_row_slice(3, 2, &[1., 0., 0., 1., 1., 1.]);
    assert!(matches!(
        prepare_aeq_and_beq(a.clone(), DVector::from_vec(vec![1., 2., 3.])),
        EqualitySystem::Unique(..)
    ));
    assert!(matches!(
        prepare_aeq_and_beq(a, DVector::from_vec(vec![1., 2., 4.])),
        EqualitySystem::Inconsistent
    ));
    assert!(matches!(
        prepare_aeq_and_beq(
            DMatrix::from_row_slice(1, 2, &[1., 1.]),
            DVector::from_vec(vec![1.])
        ),
        EqualitySystem::Multiple(..)
    ));
}
#[test]
fn square_exhaustive() {
    let n = vec![vec![1., 0.], vec![-1., 0.], vec![0., 1.], vec![0., -1.]];
    let r = ehz_from_facets(&n, &[1.; 4]).unwrap();
    assert!((r.capacity - 4.).abs() < 1e-10);
    assert_eq!(r.supports, 11);
    let score: f64 = (0..r.order.len())
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .map(|(i, j)| r.beta[i] * r.beta[j] * omega(&n[r.order[i]], &n[r.order[j]]))
        .sum();
    assert!((score - r.minimum).abs() < 1e-12);
}

#[test]
fn homogeneity_and_row_scaling() {
    let n = vec![vec![1., 0.], vec![-1., 0.], vec![0., 1.], vec![0., -1.]];
    for scale in [1e-9, 1e-6, 1., 1e6, 1e9] {
        let r = ehz_from_facets(&n, &[scale; 4]).unwrap();
        assert!((r.capacity / (4. * scale * scale) - 1.).abs() < 1e-9);
        assert!((r.maximum + r.minimum).abs() < 1e-9 / scale / scale);
        assert_eq!(
            r.maximizing_order,
            r.order.iter().rev().copied().collect::<Vec<_>>()
        );
        let scaled: Vertices = n
            .iter()
            .map(|r| r.iter().map(|x| x * scale).collect())
            .collect();
        assert!((ehz_from_facets(&scaled, &[scale; 4]).unwrap().capacity - 4.).abs() < 1e-9);
    }
}
#[test]
fn rejects_unbounded_and_zero_normals() {
    let n = vec![
        vec![1., 0., 0., 0.],
        vec![-1., 0., 0., 0.],
        vec![0., 0., 1., 0.],
        vec![0., 0., -1., 0.],
        vec![0., 1., 0., 0.],
        vec![0., 0., 0., 1.],
    ];
    assert!(ehz_from_facets(&n, &[1.; 6])
        .unwrap_err()
        .contains("unbounded"));
    assert!(ehz_from_facets(&vec![vec![1., 0.], vec![-1., 0.], vec![0., 0.]], &[1.; 3]).is_err());
}
