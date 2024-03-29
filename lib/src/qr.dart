library grizzly.linalg.qr;

import 'dart:math';
import 'package:grizzly/grizzly.dart';

QR qr(Num2D matrix) => QR.compute(matrix);

/// The reduced QR-decomposition of a matrix.
///
/// Decomposes an M x N matrix `A`, with `M >= N`, into an M x N orthogonal
/// matrix `Q` and an N x N upper rectangular matrix `R`, such that
/// `A = QR`.
///
/// The primary use of the reduced QR-decomposition is in the least squares
/// solution of non-square systems of simultaneous linear equations. This will
/// fail if the matrix is rank deficient.
class QR {
  /// QR decomposition values.
  final Double2D qr;

  /// The values on the diagonal of the upper rectangular factor.
  final Double1D rDiag;

  final Double2D q;

  final Double2D r;

  QR(this.qr, this.rDiag)
      : q = QR.orthogonalFactor(qr),
        r = QR.upperTriangularFactor(qr, rDiag);

  /// Whether or not the decomposed [Matrix] is full rank.
  bool get isFullRank {
    for (int j = 0; j < qr.numCols; j++) {
      if (rDiag[j].abs() < 0.00001) {
        return false;
      }
    }

    return true;
  }

  Double2D solve(Double2D b) {
    if (b.numRows != qr.numRows)
      throw ArgumentError('Matrix row dimensions must agree.');

    if (!isFullRank) throw UnsupportedError('Matrix is rank deficient.');

    // Copy right hand side
    final int xCols = b.numCols;
    final Double2D x = b.clone();

    // Compute Y = transpose(Q)*B
    for (int k = 0; k < qr.numCols; k++) {
      for (int j = 0; j < xCols; j++) {
        double s = 0.0;

        for (int i = k; i < qr.numRows; i++) {
          s += qr[i][k] * x[i][j];
        }

        s = -s / qr[k][k];

        for (int i = k; i < qr.numRows; i++) {
          x[i][j] += s * qr[i][k];
        }
      }
    }

    // Solve R*X = Y;
    for (int k = qr.numCols - 1; k >= 0; k--) {
      for (int j = 0; j < xCols; j++) {
        x[k][j] /= rDiag[k];
      }

      for (int i = 0; i < k; i++) {
        for (int j = 0; j < xCols; j++) {
          x[i][j] -= x[k][j] * qr[i][k];
        }
      }
    }

    return x.cut(Index2D(0, 0), Index2D(qr.numCols, xCols));
  }

  /// Creates a new [ReducedQRDecomposition] for the [matrix].
  factory QR.compute(Num2D matrix) {
    final qr = matrix.toDouble();
    final rDiag = Double1D.filled(matrix.numCols, 0);

    final int numRows = matrix.numRows;
    final int numCols = matrix.numCols;

    // Main loop.
    for (int k = 0; k < numCols; k++) {
      // Compute 2-norm of k-th column
      double nrm = 0.0;

      for (int r = k; r < numRows; r++) {
        nrm = sqrt(pow(nrm, 2) + pow(qr[r][k], 2));
      }

      if (nrm != 0.0) {
        // Form k-th Householder vector.
        if (qr[k][k] < 0) {
          nrm = -nrm;
        }

        for (int r = k; r < numRows; r++) {
          qr[r][k] /= nrm;
        }

        qr[k][k] += 1.0;

        // Apply transformation to remaining columns.
        for (int j = k + 1; j < numCols; j++) {
          double s = 0.0;

          for (int i = k; i < numRows; i++) {
            s += qr[i][k] * qr[i][j];
          }

          s = -s / qr[k][k];

          for (int i = k; i < numRows; i++) {
            qr[i][j] += s * qr[i][k];
          }
        }
      }

      rDiag[k] = -nrm;
    }

    return QR(qr, rDiag);
  }

  /// This [ReducedQRDecomposition]'s Householder matrix.
  ///
  /// Lower trapezoidal [Matrix] whose columns define the reflections.
  static Double2D householderMatrix(Double2D qr) {
    final values = qr.shaped(0.0);

    for (int i = 0; i < qr.numRows; i++) {
      for (int j = 0; j < qr.numCols; j++) {
        if (i >= j) {
          values[i][j] = qr[i][j];
        } else {
          values[i][j] = 0.0;
        }
      }
    }

    return values;
  }

  /// Computes and returns the upper triangular factor R.
  static Double2D upperTriangularFactor(Double2D qr, Double1D rDiag) {
    final int numCols = qr.numCols;
    final values = Double.sized(qr.numCols, qr.numCols);

    for (int i = 0; i < numCols; i++) {
      for (int j = 0; j < numCols; j++) {
        if (i < j) {
          values[i][j] = qr[i][j];
        } else if (i == j) {
          values[i][j] = rDiag[i];
        } else {
          values[i][j] = 0.0;
        }
      }
    }

    return values;
  }

  /// Computes and returns the orthogonal factor Q.
  static Double2D orthogonalFactor(Double2D qr) {
    final values = qr.shaped(0.0);

    for (int k = qr.numCols - 1; k >= 0; k--) {
      for (int i = 0; i < qr.numRows; i++) {
        values[i][k] = 0.0;
      }

      if (k < qr.numRows) {
        values[k][k] = 1.0;
      }

      for (int j = k; j < qr.numCols; j++) {
        if (k < qr.numRows && qr[k][k] != 0.0) {
          double s = 0.0;

          for (int i = k; i < qr.numRows; i++) {
            s += qr[i][k] * values[i][j];
          }

          s = -s / qr[k][k];

          for (int i = k; i < qr.numRows; i++) {
            values[i][j] += s * qr[i][k];
          }
        }
      }
    }

    return values;
  }
}
