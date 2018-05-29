library grizzly.linalg.svd;

import 'dart:math' as math;
import 'package:grizzly_array/grizzly_array.dart';

/// Compute the (Moore-Penrose) pseudo-inverse of the matrix [a].
///
/// Calculates the generalized inverse of a matrix using its singular-value
/// decomposition (SVD) and including all large singular values.
///
/// [rcond] specifies cutoff for small singular values. Singular values smaller
/// (in modulus) than rcond * largest_singular_value (again, in modulus) are
/// set to zero. Broadcasts against the stack of matrices.
Double2D pinv(Numeric2DView a,
    {double rcond = 1e-15, ArrayFix<double> singularValues}) {
  SVD asvd = svd(a);
  if (singularValues != null) {
    if (singularValues.length == asvd.s.length ||
        singularValues is Array<double>) {
      singularValues.assign = asvd.s;
    } else {
      throw new ArgumentError('Cannot assign singularValues!');
    }
  }
  return asvd.pinv(rcond: rcond);
}

/// Contains the [u], [s] and [v] of a matrix.
class SVD {
  final Double2DView u;

  final Double1DView s;

  final Double2DView v;

  SVD(this.u, this.s, this.v);

  /// Compute the (Moore-Penrose) pseudo-inverse of the matrix.
  ///
  /// [rcond] specifies cutoff for small singular values. Singular values smaller
  /// (in modulus) than rcond * largest_singular_value (again, in modulus) are
  /// set to zero. Broadcasts against the stack of matrices.
  Double2D pinv({double rcond = 1e-15}) {
    final Double1D sPinv = new Double1D.shapedLike(s);
    double cutoff = rcond * s.max;
    for (int i = 0; i < sPinv.length; i++) {
      if (s[i] > cutoff)
        sPinv[i] = 1 / s[i];
      else
        sPinv[i] = 0.0;
    }
    return v.matmulDiag(sPinv).matmul(u.transpose);
  }

  Double2D solve(Double2D b, {double rcond = 1e-15}) {
    if (u.numRows != b.numRows) throw new Exception('Invalid size!');
    return pinv(rcond: rcond).matmul(b);
  }

  /// Computes and returns the original matrix [a] from [u], [s] and [v].
  Double2D get original => u.matmulDiag(s).matmul(v.transpose);

  int rank({double rcond = 1e-15}) {
    final double cutoff = rcond * s.max;
    int ret = 0;
    for (int i = 0; i < s.length; i++) {
      if (s[i] > cutoff) ret++;
    }
    return ret;
  }
}

/// Computes Singular Value Decomposition of [a]
///
/// [a] is factorized into [u], [s] and [v] such that
///
/// `a = u * s * v^T`
///
/// Reference: http://cacs.usc.edu/education/phys516/src/TB/svdcmp.c
SVD svd(Numeric2DView a) {
  final int m = a.numRows;
  final int n = a.numCols;

  if (m < n) throw new ArgumentError.value(m, 'm', 'Must be >= n!');

  Double2D u = new Double2D.fromNums(a);
  final s = new Double1D.sized(n);
  final v = new Double2D.sized(n, n);

  final rv1 = new Double1D.sized(n);

  int flag, i, j, jj, k, l, nm;
  double anorm, c, f, g, h, ts, scale, x, y, z;

  g = scale = anorm = 0.0;
  // Householder reduction to bidiagonal form
  for (i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = ts = scale = 0.0;
    if (i < m) {
      for (k = i; k < m; k++) scale += u[k][i].abs();
      if (scale != 0) {
        for (k = i; k < m; k++) {
          u[k][i] /= scale;
          ts += u[k][i] * u[k][i];
        }
        f = u[i][i];
        g = -_copySign(math.sqrt(ts), f);
        h = f * g - ts;
        u[i][i] = f - g;
        for (j = l; j < n; j++) {
          ts = 0.0;
          for (k = i; k < m; k++) ts += u[k][i] * u[k][j];
          f = ts / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
        for (k = i; k < m; k++) u[k][i] *= scale;
      }
    }
    s[i] = scale * g;
    g = ts = scale = 0.0;
    if (i < m && i != n) {
      for (k = l; k < n; k++) scale += u[i][k].abs();
      if (scale != 0) {
        for (k = l; k < n; k++) {
          u[i][k] /= scale;
          ts += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = -_copySign(math.sqrt(ts), f);
        h = f * g - ts;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++) {
          ts = 0.0;
          for (k = l; k < n; k++) ts += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += ts * rv1[k];
        }
        for (k = l; k < n; k++) u[i][k] *= scale;
      }
    }
    anorm = math.max(anorm, (s[i].abs() + rv1[i].abs()));
  }

  // Accumulation of right-hand transformations
  for (i = n - 1; i >= 0; i--) {
    if (i < n) {
      if (g != 0) {
        // Double division to avoid possible underflow
        for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        for (j = l; j < n; j++) {
          ts = 0.0;
          for (k = l; k < n; k++) ts += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += ts * v[k][i];
        }
      }
      for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  // Accumulation of left-hand transformations
  for (i = math.min(m, n) - 1; i >= 0; i--) {
    l = i + 1;
    g = s[i];
    for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g != 0) {
      g = 1.0 / g;
      for (j = l; j < n; j++) {
        ts = 0.0;
        for (k = l; k < m; k++) ts += u[k][i] * u[k][j];
        f = (ts / u[i][i]) * g;
        for (k = i; k < m; k++) u[k][j] += f * u[k][i];
      }
      for (j = i; j < m; j++) u[j][i] *= g;
    } else
      for (j = i; j < m; j++) u[j][i] = 0.0;
    ++u[i][i];
  }

  for (k = n - 1; k >= 0; k--) {
    /* Diagonalization of the bidiagonal form. */
    for (int its = 1; its <= _numIterations; its++) {
      flag = 1;
      for (l = k; l >= 0; l--) {
        /* Test for splitting. */
        nm = l - 1; /* Note that rv1[1] is always zero. */
        if ((rv1[l].abs() + anorm) == anorm) {
          flag = 0;
          break;
        }
        if (nm >= 0 && (s[nm].abs() + anorm) == anorm) break;
      }
      if (flag != 0) {
        c = 0.0; /* Cancellation of rv1[l], if l > 1. */
        ts = 1.0;
        for (i = l; i <= k; i++) {
          f = ts * rv1[i];
          rv1[i] = c * rv1[i];
          if ((f.abs() + anorm) == anorm) break;
          g = s[i];
          h = _pythag(f, g);
          s[i] = h;
          h = 1.0 / h;
          c = g * h;
          ts = -f * h;
          for (j = 0; j < m; j++) {
            y = nm >= 0 ? u[j][nm] : 0.0;
            z = u[j][i];
            if (nm >= 0) u[j][nm] = y * c + z * ts;
            u[j][i] = z * c - y * ts;
          }
        }
      }
      z = s[k];
      if (l == k) {
        /* Convergence. */
        if (z < 0.0) {
          /* Singular value is made nonnegative. */
          s[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == _numIterations)
        throw new Exception("no convergence in 30 svdcmp iterations");
      x = s[l]; /* Shift from bottom 2-by-2 minor. */
      nm = k - 1;
      y = nm >= 0 ? s[nm] : 0.0;
      g = nm >= 0 ? rv1[nm] : 0.0;
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = _pythag(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + _copySign(g, f))) - h)) / x;
      c = ts = 1.0; /* Next QR transformation: */
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = s[i];
        h = ts * g;
        g = c * g;
        z = _pythag(f, h);
        if (j >= 0) rv1[j] = z;
        c = f / z;
        ts = h / z;
        f = x * c + g * ts;
        g = g * c - x * ts;
        h = y * ts;
        y *= c;
        for (jj = 0; jj < n; jj++) {
          if (j >= 0) x = v[jj][j];
          z = v[jj][i];
          if (j >= 0) v[jj][j] = x * c + z * ts;
          v[jj][i] = z * c - x * ts;
        }
        z = _pythag(f, h);
        if (j >= 0) s[j] = z; // Rotation can be arbitrary if z = 0
        if (z != 0) {
          z = 1.0 / z;
          c = f * z;
          ts = h * z;
        }
        f = c * g + ts * y;
        x = c * y - ts * g;
        for (jj = 0; jj < m; jj++) {
          if (j >= 0) y = u[jj][j];
          z = u[jj][i];
          if (j >= 0) u[jj][j] = y * c + z * ts;
          u[jj][i] = z * c - y * ts;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      s[k] = x;
    }
  }

  return new SVD(u, s, v);
}

double _copySign(double a, double b) => b >= 0.0 ? a.abs() : -a.abs();

double _pythag(double a, double b) {
  double at = a.abs(), bt = b.abs(), ct, result;

  if (at > bt) {
    ct = bt / at;
    result = at * math.sqrt(1.0 + ct * ct);
  } else if (bt > 0.0) {
    ct = at / bt;
    result = bt * math.sqrt(1.0 + ct * ct);
  } else
    result = 0.0;
  return (result);
}

const int _numIterations = 3000;
