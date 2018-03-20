/// Reference: http://cacs.usc.edu/education/phys516/src/TB/svdcmp.c
library grizzly.linalg.svd;

import 'dart:math' as math;
import 'package:grizzly_array/grizzly_array.dart';

class SVD {
  final Double2D u;

  final Double1D s;

  final Double2D v;

  SVD(Numeric2DView a)
      : u = new Double2D.fromNum(a),
        s = new Double1D.sized(a.numCols),
        v = new Double2D.shapedLike(a) {
    final int m = a.numRows;
    final int n = a.numCols;

    if (m < n) throw new ArgumentError.value(m, 'm', 'Must be >= n!');

    final rv1 = new Double1D.sized(n);

    int flag, i, its, j, jj, k, l, nm;
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
          g = -copySign(math.sqrt(ts), f);
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
          g = -copySign(math.sqrt(ts), f);
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
      for (its = 1; its <= numIterations; its++) {
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
        if (its == numIterations)
          throw new Exception("no convergence in 30 svdcmp iterations");
        x = s[l]; /* Shift from bottom 2-by-2 minor. */
        nm = k - 1;
        y = nm >= 0 ? s[nm] : 0.0;
        g = nm >= 0 ? rv1[nm] : 0.0;
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = _pythag(f, 1.0);
        f = ((x - z) * (x + z) + h * ((y / (f + copySign(g, f))) - h)) / x;
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
  }

  static double copySign(double a, double b) => b >= 0.0 ? a.abs() : -a.abs();

  static double _pythag(double a, double b) {
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

  static const int numIterations = 3000;
}

/*
  SVD1(Numeric2DView a)
      : u = new Double2D.fromNum(a),
        s = new Double1D.sized(a.numCols + 1),
        v = new Double2D.shapedLike(a) {
    final int m = a.numRows;
    final int n = a.numCols;

    u.row.insertScalar(0, 0.0);
    u.col.insertScalar(0, 0.0);

    v.row.insertScalar(0, 0.0);
    v.col.insertScalar(0, 0.0);

    if (m < n) throw new ArgumentError.value(m, 'm', 'Must be >= n!');

    final rv1 = new Double1D.sized(n + 1);

    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, ts, scale, x, y, z;

    g = scale = anorm = 0.0;
    // Householder reduction to bidiagonal form
    for (i = 1; i <= n; i++) {
      l = i + 1;
      rv1[i] = scale * g;
      g = ts = scale = 0.0;
      if (i <= m) {
        for (k = i; k <= m; k++) scale += u[k][i].abs();
        if (scale != 0) {
          for (k = i; k <= m; k++) {
            u[k][i] /= scale;
            ts += u[k][i] * u[k][i];
          }
          f = u[i][i];
          g = -copySign(math.sqrt(ts), f);
          h = f * g - ts;
          u[i][i] = f - g;
          for (j = l; j <= n; j++) {
            ts = 0.0;
            for (k = i; k <= m; k++) ts += u[k][i] * u[k][j];
            f = ts / h;
            for (k = i; k <= m; k++) u[k][j] += f * u[k][i];
          }
          for (k = i; k <= m; k++) u[k][i] *= scale;
        }
      }
      s[i] = scale * g;
      g = ts = scale = 0.0;
      if (i <= m && i != n) {
        for (k = l; k <= n; k++) scale += u[i][k].abs();
        if (scale != 0) {
          for (k = l; k <= n; k++) {
            u[i][k] /= scale;
            ts += u[i][k] * u[i][k];
          }
          f = u[i][l];
          g = -copySign(math.sqrt(ts), f);
          h = f * g - ts;
          u[i][l] = f - g;
          for (k = l; k <= n; k++) rv1[k] = u[i][k] / h;
          for (j = l; j <= m; j++) {
            ts = 0.0;
            for (k = l; k <= n; k++) ts += u[j][k] * u[i][k];
            for (k = l; k <= n; k++) u[j][k] += ts * rv1[k];
          }
          for (k = l; k <= n; k++) u[i][k] *= scale;
        }
      }
      anorm = math.max(anorm, (s[i].abs() + rv1[i].abs()));
    }

    for (i = n; i >= 1; i--) {
      /* Accumulation of right-hand transformations. */
      if (i < n) {
        if (g != 0) {
          for (j = l;
              j <= n;
              j++) /* Double division to avoid possible underflow. */
            v[j][i] = (u[i][j] / u[i][l]) / g;
          for (j = l; j <= n; j++) {
            ts = 0.0;
            for (k = l; k <= n; k++) ts += u[i][k] * v[k][j];
            for (k = l; k <= n; k++) v[k][j] += ts * v[k][i];
          }
        }
        for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
      }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }
    for (i = math.min(m, n); i >= 1; i--) {
      /* Accumulation of left-hand transformations. */
      l = i + 1;
      g = s[i];
      for (j = l; j <= n; j++) u[i][j] = 0.0;
      if (g != 0) {
        g = 1.0 / g;
        for (j = l; j <= n; j++) {
          ts = 0.0;
          for (k = l; k <= m; k++) ts += u[k][i] * u[k][j];
          f = (ts / u[i][i]) * g;
          for (k = i; k <= m; k++) u[k][j] += f * u[k][i];
        }
        for (j = i; j <= m; j++) u[j][i] *= g;
      } else
        for (j = i; j <= m; j++) u[j][i] = 0.0;
      ++u[i][i];
    }

    for (k = n; k >= 1; k--) {
      /* Diagonalization of the bidiagonal form. */
      for (its = 1; its <= 30; its++) {
        flag = 1;
        for (l = k; l >= 1; l--) {
          /* Test for splitting. */
          nm = l - 1; /* Note that rv1[1] is always zero. */
          if ((rv1[l].abs() + anorm) == anorm) {
            flag = 0;
            break;
          }
          if ((s[nm].abs() + anorm) == anorm) break;
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
            for (j = 1; j <= m; j++) {
              y = u[j][nm];
              z = u[j][i];
              u[j][nm] = y * c + z * ts;
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
            for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
          }
          break;
        }
        if (its == 30)
          throw new Exception("no convergence in 30 svdcmp iterations");
        x = s[l]; /* Shift from bottom 2-by-2 minor. */
        nm = k - 1;
        y = s[nm];
        g = rv1[nm];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = _pythag(f, 1.0);
        f = ((x - z) * (x + z) + h * ((y / (f + copySign(g, f))) - h)) / x;
        c = ts = 1.0; /* Next QR transformation: */
        for (j = l; j <= nm; j++) {
          i = j + 1;
          g = rv1[i];
          y = s[i];
          h = ts * g;
          g = c * g;
          z = _pythag(f, h);
          rv1[j] = z;
          c = f / z;
          ts = h / z;
          f = x * c + g * ts;
          g = g * c - x * ts;
          h = y * ts;
          y *= c;
          for (jj = 1; jj <= n; jj++) {
            x = v[jj][j];
            z = v[jj][i];
            v[jj][j] = x * c + z * ts;
            v[jj][i] = z * c - x * ts;
          }
          z = _pythag(f, h);
          s[j] = z; /* Rotation can be arbitrary if z = 0. */
          if (z != 0) {
            z = 1.0 / z;
            c = f * z;
            ts = h * z;
          }
          f = c * g + ts * y;
          x = c * y - ts * g;
          for (jj = 1; jj <= m; jj++) {
            y = u[jj][j];
            z = u[jj][i];
            u[jj][j] = y * c + z * ts;
            u[jj][i] = z * c - y * ts;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        s[k] = x;
      }
    }
  }
*/
