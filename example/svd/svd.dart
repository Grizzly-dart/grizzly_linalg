import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

Double2D flipv(Double2D other) {
  Double2D ret = new Double2D.shapedLike(other);
  for (int i = 0; i < other.numRows; i++) {
    for (int j = 0; j < other.numCols; j++) {
      ret[i][j] = other[i][other.numCols - j - 1];
    }
  }
  return ret;
}

Double2D flip(Double2D other) {
  Double2D ret = new Double2D.shapedLike(other);
  for (int i = 0; i < other.numRows; i++) {
    for (int j = 0; j < other.numCols; j++) {
      ret[i][j] = other[other.numRows - i - 1][other.numCols - j - 1];
    }
  }
  return ret;
}

main() {
  Double2D x = new Int2D([
    [1, 1, 5],
    [1, 2, 8],
    [1, 3, 4],
    [1, 4, 10],
    [1, 5, 11],
  ]).toDouble();

  /*
  x = new Int2D([
    [1, 1, 0, 0, 0],
    [1, 2, 0, 0, 0],
    [1, 3, 0, 0, 0],
    [1, 4, 0, 0, 0],
    [1, 5, 0, 0, 0],
  ]).toDouble();
*/

  print('x');
  print(x);

  Double2D y = new Double2D.aCol(x.col[0] + (x.col[1] * 5));
  // print('y');
  // print(y);

  final xqr = svd(x);

  print(xqr.u.shape);
  print(xqr.s.shape);
  print(xqr.v.shape);

  Double2D u = xqr.u;
  Double2D s = xqr.s.diagonal();
  Double2D v = xqr.v;

  v = v.transpose;

  print(u);
  print(s);
  print(v);

  print(u.matmul(s).matmul(v));

  /*
  print('X: q r');
  print(xqr.q);
  print(xqr.r);

  print('');
  print(xqr.q * xqr.r);
  */

  /*
  final y2d = new Double2D.(y);
  print(y.shape);
  print(y2d.shape);
  print(xqr.solve(y2d));

  print(xqr.solve(y));
  */
}
