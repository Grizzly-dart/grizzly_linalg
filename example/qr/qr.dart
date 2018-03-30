import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

/// shape(a) = 5 x 2
/// shape(aT) = 2 x 5
/// shape(aT * a) = 2 x 2
///
/// shape(x) = 2 x 1
///
/// shape(b) = 5 x 1
/// shape(aT * b) = 2 x 1

main() {
  final x = array2D([
    [1, 1],
    [1, 2],
    [1, 3],
    [1, 4],
    [1, 5],
  ]);

  print('x');
  print(x);

  Double2D y = new Double2D.aCol(x.col[0] + (x.col[1] * 5));
  print('y');
  print(y);

  final xqr = qr(x);
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
  */

  print(xqr.solve(y));
}
