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

  final y = x.col[0] * 5;
  print(y);

  final xqr = qr(x);
  print(xqr.q);
  print(xqr.r);

  print(xqr.q.shape);
  print(xqr.r.shape);

  print(xqr.q * xqr.r);

  /*
  final y2d = new Double2D.columns(y);
  print(y.shape);
  print(y2d.shape);
  print(xqr.solve(y2d));
  */
}