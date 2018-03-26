import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  final x = array2D([
    [1],
    [2],
    [3],
    [4],
    [5],
  ]);

  print('x');
  print(x);

  final y = (x * 5) /* TODO  + [1, 1, 1, 1, 1] */;
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