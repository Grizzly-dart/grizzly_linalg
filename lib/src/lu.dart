library grizzly.linalg.lu;

import 'package:grizzly/grizzly.dart';

LU lu(Double2D matrix) => LU.compute(matrix);

Double2D solve(Num2D a, Num2D b) => LU.compute(a).solve(b);

double det(Num2D a) => LU.compute(a).determinant;

/// The lower-upper factor decomposition of a [Matrix], with partial pivoting.
///
/// Lower-upper factor decomposition with partial pivoting of an M x N matrix
/// `A`, results in 3 matrices:
///
/// - `L`: the lower factor matrix. An M x N [Matrix] with all zero's above the
///   diagonal.
/// - `U`: the upper factor matrix. An N x N [Matrix] with all zero's below the
///   diagonal.
/// - `P`: the pivot matrix. An M x M permutation [Matrix].
///
/// Such that `PA = LU`.
///
/// The primary use of the lower-upper decomposition is in the solution of
/// square systems of simultaneous linear equations. This will fail if the
/// [Matrix] is non-singular. Pivoting reduces the impact of rounding errors.
class LU {
  /// LU decomposition values.
  Double2D lu;

  /// The pivot vector.
  ///
  /// Used to keep track of where to place 1's in the pivot matrix.
  ///
  /// Each row in the pivot matrix consists of zero's and one 1. The values
  /// in the pivot vector track in which column to place the one (zero indexed).
  /// A pivot vector of `(0, 2, 3, 1)` thus corresponds to the following pivot
  /// matrix:
  ///
  ///     1 0 0 0
  ///     0 0 1 0
  ///     0 0 0 1
  ///     0 1 0 0
  ///
  Int1D piv;

  /// The pivot sign.
  num pivotSign;

  LU(this.lu, this.piv, this.pivotSign);

  /// Creates a new [LU] for the [matrix].
  factory LU.compute(Num2DView matrix) {
    final lu = matrix.toDouble();
    final int numRows = lu.numRows;
    final int numCols = lu.numCols;
    int pivotSign = 1;
    final piv = List<int>.generate(numRows, (i) => i);

    // Outer loop.
    for (int j = 0; j < numCols; j++) {
      // Find pivot
      int p = j;

      for (int i = j + 1; i < numRows; i++) {
        if (lu[i][j].abs() > lu[p][j].abs()) {
          p = i;
        }
      }

      // Exchange pivot if necessary
      if (p != j) {
        for (int k = 0; k < numCols; k++) {
          final double t = lu[p][k];
          lu[p][k] = lu[j][k];
          lu[j][k] = t;
        }

        final int k = piv[p];
        piv[p] = piv[j];
        piv[j] = k;

        pivotSign = -pivotSign;
      }

      // Compute multipliers.
      if (j < numRows && lu[j][j] != 0.0) {
        for (int i = j + 1; i < numRows; i++) {
          lu[i][j] /= lu[j][j];

          for (int k = j + 1; k < numCols; k++) {
            lu[i][k] -= lu[i][j] * lu[j][k];
          }
        }
      }
    }

    return LU(lu, piv, pivotSign);
  }

  /// Whether is not the decomposed [Matrix] is non-singular.
  ///
  /// A non-singular [Matrix] has an inverse and a non-zero determinant.
  ///
  /// Throws an [UnsupportedError] if the decomposed [Matrix] is not square.
  bool get isSingular {
    if (!lu.isSquare) throw UnsupportedError('Matrix is not square.');

    for (int j = 0; j < lu.numCols; j++) {
      if (lu[j][j] == 0.0) return true;
    }

    return false;
  }

  /// This [PivotingLUDecomposition]'s lower factor.
  ///
  /// A [Matrix] with all zero's above the diagonal.
  Double2D get l {
    final values = lu.shaped<double>(0);

    for (int i = 0; i < lu.numRows; i++) {
      for (int j = 0; j < lu.numCols; j++) {
        if (i > j) {
          values[i][j] = lu[i][j];
        } else if (i == j) {
          values[i][j] = 1.0;
        }
      }
    }

    return values;
  }

  /// This [PivotingLUDecomposition]'s upper factor.
  ///
  /// A [Matrix] with all zero's below the diagonal.
  Double2D get u {
    final values = Double.filled2D(lu.numCols, lu.numCols);

    for (int i = 0; i < lu.numCols; i++) {
      for (int j = 0; j < lu.numCols; j++) {
        if (i <= j) {
          values[i][j] = lu[i][j];
        }
      }
    }

    return values;
  }

  /// This [PivotingLUDecomposition]'s pivot matrix.
  ///
  /// A permutation matrix.
  Double2D get p {
    final values = Double.filled2D(lu.numRows, lu.numRows);

    for (int i = 0; i < lu.numRows; i++) {
      for (int j = 0; j < lu.numRows; j++) {
        if (i == piv[j]) {
          values[i][j] = 1.0;
        }
      }
    }

    return values;
  }

  /// The decomposed [Matrix]'s determinant.
  ///
  /// Throws an [UnsupportedError] if the decomposed [Matrix] is not square.
  double get determinant {
    if (!lu.isSquare) throw UnsupportedError('Matrix must be square.');

    double determinant = pivotSign.toDouble();

    for (int j = 0; j < lu.numCols; j++) determinant *= lu[j][j];

    return determinant;
  }

  /// Solves `AX=B` for X, where `A` is the decomposed [Matrix] and [b] the
  /// given [Matrix].
  ///
  /// Throws an [ArgumentError] if the row dimensions of `A` and [b] do not
  /// match.
  ///
  /// Throws an [UnsupportedError] if `A` is not square.
  ///
  /// Throws an [UnsupportedError] if `A` is singular (not invertible).
  Double2D solve(Num2D b) {
    if (b.numRows != lu.numRows)
      throw ArgumentError('Matrix row dimensions must agree.');

    if (isSingular) throw UnsupportedError('Matrix is singular.');

    final xCols = b.numCols;
    final x = Double.filled2D(lu.numCols, xCols);

    // Copy right hand side with pivoting
    {
      int c = 0;
      for (int row in piv) {
        x[c++] = b[row].toDouble();
      }
    }

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < lu.numCols; k++) {
      for (int i = k + 1; i < lu.numCols; i++) {
        for (int j = 0; j < xCols; j++) {
          x[i][j] -= x[k][j] * lu[i][k];
        }
      }
    }

    // Solve U*X = Y;
    for (int k = lu.numCols - 1; k >= 0; k--) {
      for (int j = 0; j < xCols; j++) {
        x[k][j] /= lu[k][k];
      }

      for (int i = 0; i < k; i++) {
        for (int j = 0; j < xCols; j++) {
          x[i][j] -= x[k][j] * lu[i][k];
        }
      }
    }

    return x;
  }
}
