import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  final a = new Double2D([
    [1.0, 2.0],
    [2.0, 4.0],
    [3.0, 6.0]
  ]);

  print('A');
  print(a);

  final aqr = qr(a);

  print('Q');
  print(aqr.q);
  print('R');
  print(aqr.r);

  print('Full rank? ${aqr.isFullRank}');

  /*
  print('A = Q * R');
  print(aqr.q.matmul(aqr.r));

  Double2D b = new Double2D.aCol(a.col[0] + (a.col[1] * 5));
  print('B');
  print(b);

  print('x = A^-1 * b');
  print(aqr.solve(b));
  */
}