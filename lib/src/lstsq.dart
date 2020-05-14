library grizzly.linalg.lstsq;

import 'package:grizzly_array/grizzly_array.dart';

/// Return the least-squares solution to a linear matrix equation
///
/// Solves the equation a x = b by computing a vector x that minimizes the
/// Euclidean 2-norm || b - a x ||^2. The equation may be under-, well-, or
/// over- determined (i.e., the number of linearly independent rows of a can be
/// less than, equal to, or greater than its number of linearly independent
/// columns).
Double1D lstsqBGD(Numeric2DView x, Numeric1DView y,
        {double learningRate: 1e-4,
        int maxIterations: 200,
        Iterable<double> initParams}) =>
    (new BatchLeastSquareGradientDescent(
      x,
      y,
      learningRate: learningRate,
      maxIterations: maxIterations,
      initParams: initParams,
    )..learn())
        .params;

Double1D lstsqSGD(Numeric2DView x, Numeric1DView y,
        {double learningRate: 1e-4,
        int maxIterations: 200,
        Iterable<double> initParams}) =>
    (new StochasticLeastSquareGradientDescent(
      x,
      y,
      learningRate: learningRate,
      maxIterations: maxIterations,
      initParams: initParams,
    )..learn())
        .params;

abstract class LeastSquareGradientDescent {
  /// Learning rate used in gradient descent
  double get learningRate;

  /// Maximum iterations
  int get maxIterations;

  /// (N, 1) parameters being estimated
  Double1DFix get params;

  /// (M, N) exogenous or independent matrix
  Numeric2DView get x;

  /// (M, 1) endogenous of dependent matrix
  Numeric1DView get y;

  /// Finds the value of the hypothesis function given the current parameter
  /// vector θ ([params]) and a sample x ([row])
  // TODO normalize argument?
  double predict(Iterable<num> row) => params.dot(row);

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
    return sumSquaredError / (2 * x.numRows);
  }
}

class BatchLeastSquareGradientDescent extends LeastSquareGradientDescent {
  /// Learning rate used in gradient descent
  final double learningRate;

  /// Maximum iterations
  final int maxIterations;

  /// (N, 1) parameters being estimated
  final Double1DFix params;

  /// (M, N) exogenous or independent matrix
  final Numeric2DView x;

  /// (M, 1) endogenous of dependent matrix
  final Numeric1DView y;

  BatchLeastSquareGradientDescent(this.x, this.y,
      {this.learningRate: 1e-4,
      this.maxIterations: 200,
      Iterable<double> initParams})
      : params = initParams == null
            ? new Double1DFix.sized(x.numCols)
            : new Double1DFix(initParams) {
    // Validate
    if (x.numRows != y.length) {
      throw new Exception('x and y must have same number of samples!');
    }
    if (x.numCols != params.length) {
      throw new Exception('x and params must have same number of features!');
    }
  }

  /// Performs least-square estimation through batch gradient descent
  void learn() {
    final theta = new Double1D.sized(params.length);
    for (int i = 0; i < maxIterations; i++) {
      for (int j = 0; j < params.length; j++) {
        theta[j] = params[j] + learningRate * dj(j);
      }

      // Update params
      for (int j = 0; j < params.length; j++) {
        final double newThetaJ = theta[j];

        if (newThetaJ.isInfinite || newThetaJ.isNaN) {
          throw new Exception('Learning diverged!');
        }

        params[j] = newThetaJ;
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
    return sum / x.numCols;
  }
}

class StochasticLeastSquareGradientDescent extends LeastSquareGradientDescent {
  /// Learning rate used in gradient descent
  final double learningRate;

  /// Maximum iterations
  final int maxIterations;

  final Double1DFix params;

  final Numeric2DView x;

  final Numeric1DView y;

  StochasticLeastSquareGradientDescent(this.x, this.y,
      {this.learningRate: 1e-4,
      this.maxIterations: 800,
      Iterable<double> initParams})
      : params = initParams == null
            ? new Double1DFix.sized(x.numCols)
            : new Double1DFix(initParams) {
    // Validate
    if (x.numRows != y.length) {
      throw new Exception('x and y must have same number of samples!');
    }
    if (x.numCols != params.length) {
      throw new Exception('x and params must have same number of features!');
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
    return gradient;
  }

  /// Performs least-square estimation through stochastic gradient descent
  void learn() {
    final theta = new Double1D.sized(params.length);
    for (int i = 0; i < maxIterations; i++) {
      for (int i = 0; i < x.numRows; i++) {
        for (int j = 0; j < params.length; j++) {
          theta[j] = params[j] + learningRate * dij(i, j);
        }

        // Update params
        for (int j = 0; j < params.length; j++) {
          final double newThetaJ = theta[j];

          if (newThetaJ.isInfinite || newThetaJ.isNaN) {
            throw new Exception('Learning diverged!');
          }

          params[j] = newThetaJ;
        }
      }
    }
  }
}
