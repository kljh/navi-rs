
/*

En dynamique des fluides, un écoulement est potentiel lorsque son champ des vitesses
est le gradient d'une fonction scalaire, appelée potentiel des vitesses.

Puisque le rotationnel d'un gradient est toujours égal à zéro, un écoulement potentiel est toujours irrotationnel.
Les écoulements potentiels servent le plus souvent à décrire des écoulements de fluides parfaits, c'est-à-dire des écoulements où la viscosité peut être négligée, parce qu'un écoulement irrotationnel le reste tant que la viscosité est négligeable (équation d'Euler avec l'hypothèse que le champ de forces extérieures dérive d'un potentiel).

Si l'écoulement est incompressible, la divergence de v est nulle.
Le potentiel des vitesses est alors une solution de l'équation de Laplace : laplacien(v) = 0;

*/

use std::isize;
use std::usize;

use nalgebra::DMatrix;
use nalgebra::DVector;
use nalgebra::Vector2;

fn r2(x: f64, y: f64) -> f64 {
    // are the two below equivalent ?
    // x.powi(2) + y.powi(2)
    x*x + y*y
}

pub trait Solid {
    fn is_solid(&self, x: f64, y: f64) -> bool;
    fn normal_direction(&self, x: f64, y: f64) -> Vector2<f64>;
}

pub struct FreeFlow {}

impl Solid for FreeFlow {
    fn is_solid(&self, _x: f64, _y: f64) -> bool {
        false
    }
    fn normal_direction(&self, _x: f64, _y: f64) -> Vector2<f64> {
        Vector2::new(0.0, 0.0) // not used
        // panic!("Free flow does not have a normal direction");
    }
}

pub struct Cylinder{
    pub center: Vector2<f64>,
    pub radius2: f64
}

impl Solid for Cylinder {
    fn is_solid(&self, x: f64, y: f64) -> bool {
        r2(x-self.center[0], y-self.center[1]) < self.radius2
    }
    fn normal_direction(&self, x: f64, y: f64) -> Vector2<f64> {
        let er = Vector2::new(x-self.center[0], y-self.center[1]);
        er
    }
}

/*

pub struct CompositeSolid {
    pub solids: Vec<Box<dyn Solid>>,
}

impl Solid for CompositeSolid {
    fn is_solid(&self, x: f64, y: f64) -> bool {
        self.solids.iter().any(|s| s.is_solid(x, y))
    }
    fn normal_direction(&self, x: f64, y: f64) -> Vector2<f64> {
        for s in &self.solids {
            if s.is_solid(x, y) {
                return s.normal_direction(x, y);
            }
        }
        Vector2::new(0.0, 0.0) // not used
        // panic!("Point is not solid, no normal direction");
    }

}
*/

enum LimitCondition {
    SolidPoint,
    FreeFlow,
    TangentialFlow { normal_direction: Vector2<f64> },
}

// limit condition for the flow for a cylinder
// returns None if it's an interior (or solid) point
// or Some(cond) if it's a boundary point
fn limit_condition(i: usize, j: usize, geom: &dyn Solid) -> LimitCondition {
    let x = i as f64;
    let y = j as f64;
    let dx = 1.0;
    let dy = 1.0;

    if geom.is_solid(x, y) {
        LimitCondition::SolidPoint
    }
    else if
        geom.is_solid(x+dx, y) || geom.is_solid(x-dx, y) ||
        geom.is_solid(x, y+dy) || geom.is_solid(x, y-dy)
    {
        let ndir = geom.normal_direction(x, y);
        LimitCondition::TangentialFlow { normal_direction: ndir }
    } else {
        LimitCondition::FreeFlow
    }
}

// Résolution par difference finies sur grille orthonormale nxm
pub fn finite_diff(n: usize, m: usize, geom: &dyn Solid) -> DMatrix<f64>
{
    let mut a : DMatrix<f64> = DMatrix::from_element(n*m, n*m, 0.0);
    let mut b : DVector<f64> = DVector::from_element(n*m, 0.0);

    // set zeros
    for i in 0..n {
        for j in 0..m {
            // row major
            let k = i*m + j;
            a[(k, k)] = if i==j {1.0} else {0.0};
            b[k] = 0.0;
        }
    }

    let potentiel = | _i:isize, j: isize| {
        50.+10.*(j as f64)
    };

    let mut set_weigth = | ui:usize, uj: usize, di: isize, dj: isize, w: f64, k0: usize| {
        let si = ui as isize + di;
        let sj = uj as isize + dj;

        if si<0 || si == (n as isize) || sj<0 || sj == (m as isize)  {
            let y = - w * potentiel(si, sj);
            b[k0] = b[k0] + y;
            println!("  set weight => {ui}+{di} {uj}+{dj} => {si} {sj} => b[k0={k0}]=-{w}*limit");
        } else {
            let k1 = si * m as isize + sj;
            a[(k0, k1 as usize)] = w;
            println!("  set weight => {ui}+{di} {uj}+{dj} => {si} {sj} => A[k0={k0},k1={k1}]={w}");
        }
    };

    for i in 0..n {
        for j in 0..m {
            let k0 = i*m +j;
            println!("k0 {k0}");
            match limit_condition(i, j, geom)
            {
                LimitCondition::FreeFlow => {
                    // Interior point, Laplacien is zero
                    set_weigth(i, j, 0, 0, 1.0, k0);
                    set_weigth(i, j, -1, 0, -0.25, k0);
                    set_weigth(i, j, 1, 0, -0.25, k0);
                    set_weigth(i, j, 0, -1, -0.25, k0);
                    set_weigth(i, j, 0, 1, -0.25, k0);
                },
                LimitCondition::TangentialFlow { normal_direction: ndir } => {
                    // Boundary point, flow is tangent to the boundary,
                    // so potential is constant along the normal direction to the boundary.

                    // |w[0]| + |w[1]| = 1.0, so we can split the weight in two parts for the two neighbors in the normal direction
                    let wght = ndir / (ndir[0].abs() + ndir[1].abs());

                    println!("Boundary point {}, {} ndir {} wght {}", i, j, ndir, wght);
                    set_weigth(i, j, 0, 0, 1.0, k0);
                    set_weigth(i, j, if wght[0]>0.0 {1} else {-1}, 0, -(wght[0].abs()), k0);
                    set_weigth(i, j, 0, if wght[1]>0.0 {1} else {-1}, -(wght[1].abs()), k0);
                },
                LimitCondition::SolidPoint => {
                    // Solid point, potential is zero
                    set_weigth(i, j, 0, 0, 19.87, k0);  // !! should set to medium potential for display reasons.

                    // a[(k0, k0)] = 1.0;
                    // b[k0] = 50.0;
                },
            }
        }
    }

    println!("A {}", &a);
    println!("b {}", &b);

    let dcp = a.lu();
    let x = dcp.solve(&b).unwrap();

    println!("x {:?}", &x);

    let res : DMatrix<f64> = DMatrix::from_fn(n, m, | i, j| {
        let k = i*m + j;
        x[k]
    });

    println!("res {}", &res);
    res
}


#[cfg(test)]
mod tests {
    use crate::pflow::{*};

    #[test]
    fn test_pflow() {
        let geom = FreeFlow {};
        finite_diff(4, 4, &geom);
        assert!(true);
    }

    #[test]
    fn test_pflow_cyl_5() {
        let n = 5;
        let c = ( n/2 ) as f64;
        let r2 = 0.5;
        let geom = Cylinder { center: Vector2::new(c, c), radius2: r2 };

        finite_diff(n, n, &geom);
        assert!(true);
    }

    #[test]
    fn test_pflow_cyl_9() {
        let n = 9;
        let c = ( n/2 ) as f64;
        let geom = Cylinder { center: Vector2::new(c, c), radius2: c };

        finite_diff(n, n, &geom);
        assert!(true);
    }

    #[test]
    fn test_pflow_cyl_16() {
        let n = 16;
        let c = ( n/2 ) as f64;
        let geom = Cylinder { center: Vector2::new(c, c), radius2: c };

        finite_diff(n, n, &geom);
        assert!(true);
    }
}
