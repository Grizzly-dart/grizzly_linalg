library grizzly.linalg.lstsq;

import 'package:grizzly/grizzly.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

/// Return the least-squares solution to a linear matrix equation
///
/// Solves the equation a x = b by computing a vector x that minimizes the
/// Euclidean 2-norm || b - a x ||^2. The equation may be under-, well-, or
/// over- determined (i.e., the number of linearly independent rows of a can be
/// less than, equal to, or greater than its number of linearly independent
/// columns).
List<double> lstsqBGD(List<List<num>> x, List<num> y,
        {double learningRate: 1e-4,
        int maxIterations: 200,
        List<double>? initParams}) =>
    (LeastSquareBGD(
      x,
      y,
      learningRate: learningRate,
      maxIterations: maxIterations,
      initCoeff: initParams,
    )..learn())
        .coeff;

List<double> lstsqSGD(List<List<num>> x, List<num> y,
        {double learningRate: 1e-4,
        int maxIterations: 200,
        List<double>? initParams}) =>
    (LeastSquareSGD(
      x,
      y,
      learningRate: learningRate,
      maxIterations: maxIterations,
      initCoeff: initParams,
    )..learn())
        .coeff;

abstract class LeastSquareGD {
  /// Learning rate used in gradient descent
  double get learningRate;

  /// Maximum iterations
  int get maxIterations;

  /// (N, 1) parameters being estimated
  List<double> get coeff;

  /// (M, N) exogenous or independent matrix
  List<List<num>> get x;

  /// (M, 1) endogenous of dependent matrix
  List<num> get y;

  LeastSquareRegularizer? get regularizer;

  /// Finds the value of the hypothesis function given the current parameter
  /// vector θ ([params]) and a sample x ([row])
  // TODO normalize argument?
  double predict(Iterable<num> row) => coeff.dot(row);

  /// Performs the learning
  void learn();

  /// Returns the Least Squares cost function of the given linear
  /// model. Could be useful in testing convergence.
  double costFunction() {
    double sumSquaredError = 0.0;
    for (int r = 0; r < x.numRows; r++) {
      final double prediction = predict(x[r]);
      final double error = y[r] - prediction;
      sumSquaredError += error * error;
    }
    return sumSquaredError / (2 * x.numRows) +
        (regularizer?.costFunction(coeff) ?? 0);
  }
}

class LeastSquareBGD extends LeastSquareGD {
  /// Learning rate used in gradient descent
  final double learningRate;

  /// Maximum iterations
  final int maxIterations;

  /// (N, 1) parameters being estimated
  final List<double> coeff;

  /// (M, N) exogenous or independent matrix
  final List<List<num>> x;

  /// (M, 1) endogenous of dependent matrix
  final List<num> y;

  final LeastSquareRegularizer? regularizer;

  LeastSquareBGD(this.x, this.y,
      {this.learningRate: 1e-4,
      this.maxIterations: 200,
      Num1DView? initCoeff,
      this.regularizer})
      : coeff = initCoeff == null
            ? List<double>.filled(x.numCols, 0)
            : initCoeff.toDouble() {
    // Validate
    if (x.numRows != y.length) {
      throw Exception('x and y must have same number of samples!');
    }
    if (x.numCols != coeff.length) {
      throw Exception('x and params must have same number of features!');
    }
  }

  /// Performs least-square estimation through batch gradient descent
  void learn() {
    final theta = List<double>.filled(coeff.length, 0);
    for (int i = 0; i < maxIterations; i++) {
      for (int j = 0; j < coeff.length; j++) {
        theta[j] = coeff[j] + learningRate * dj(j);
      }

      // Update params
      for (int j = 0; j < coeff.length; j++) {
        final double thetaJ = theta[j];

        if (thetaJ.isInfinite || thetaJ.isNaN) {
          throw Exception('Learning diverged!');
        }

        coeff[j] = thetaJ;
      }
    }
  }

  /// [dj] returns the partial derivative of the cost function J(θ)
  /// with respect to theta[j] where theta is the parameter vector
  /// associated with our hypothesis function Predict (upon which
  /// we are optimizing
  double dj(int j) {
    double sum = 0.0;
    for (int r = 0; r < x.numRows; r++) {
      final double prediction = predict(x[r]);
      sum += (y[r] - prediction) * x[r][j];
    }
    return sum / x.numCols + (regularizer?.derivative(coeff[j]) ?? 0);
  }
}

class LeastSquareSGD extends LeastSquareGD {
  /// Learning rate used in gradient descent
  final double learningRate;

  /// Maximum iterations
  final int maxIterations;

  final List<double> coeff;

  final List<List<num>> x;

  final List<num> y;

  final LeastSquareRegularizer? regularizer;

  LeastSquareSGD(this.x, this.y,
      {this.learningRate: 1e-4,
      this.maxIterations: 800,
      Num1DView? initCoeff,
      this.regularizer})
      : coeff = initCoeff == null
            ? List<double>.filled(x.numCols, 0)
            : initCoeff.toDouble() {
    // Validate
    if (x.numRows != y.length) {
      throw Exception('x and y must have same number of samples!');
    }
    if (x.numCols != coeff.length) {
      throw Exception('x and params must have same number of features!');
    }
  }

  /// dij returns the derivative of the cost function
  /// J(θ) with respect to the j-th parameter of
  /// the hypothesis, θ[j], for the training example
  /// x[i]. Used in Stochastic Gradient Descent.
  ///
  /// assumes that i,j is within the bounds of the
  /// data they are looking up! (because this is getting
  /// called so much, it needs to be efficient with
  /// comparisons)
  double dij(int i, int j) {
    final double prediction = predict(x[i]);
    final double gradient = (y[i] - prediction) * x[i][j];
    return gradient + (regularizer?.derivative(coeff[j]) ?? 0);
  }

  /// Performs least-square estimation through stochastic gradient descent
  void learn() {
    final theta = List<double>.filled(coeff.length, 0);
    for (int i = 0; i < maxIterations; i++) {
      for (int i = 0; i < x.numRows; i++) {
        for (int j = 0; j < coeff.length; j++) {
          theta[j] = coeff[j] + learningRate * dij(i, j);
        }

        // Update params
        for (int j = 0; j < coeff.length; j++) {
          final double newThetaJ = theta[j];

          if (newThetaJ.isInfinite || newThetaJ.isNaN) {
            throw Exception('Learning diverged!');
          }

          coeff[j] = newThetaJ;
        }
      }
    }
  }
}

abstract class LeastSquareRegularizer {
  String get name;

  double costFunction(Num1DView coefficients);

  double derivative(num coefficient);
}

class LassoRegularizer implements LeastSquareRegularizer {
  final String name = 'Lasso';

  final double lambda;

  const LassoRegularizer(this.lambda);

  double costFunction(Num1DView coefficients) =>
      lambda * coefficients.abs().sum;

  double derivative(num coefficient) => lambda * coefficient.sign;
}

class RidgeRegularizer implements LeastSquareRegularizer {
  final String name = 'Ridge';

  final double lambda;

  const RidgeRegularizer(this.lambda);

  double costFunction(Num1DView coefficients) =>
      lambda * coefficients.pow(2).sum;

  double derivative(num coefficient) => lambda * 2 * coefficient;
}
