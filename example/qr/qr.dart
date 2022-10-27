import 'package:grizzly/grizzly.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  final a = [
    [1.0, 1.0],
    [1.0, 2.0],
    [1.0, 3.0],
    [1.0, 4.0],
    [1.0, 5.0],
  ];

  print('A');
  print(a);

  final aqr = qr(a);

  print('Q');
  print(aqr.q);
  print('R');
  print(aqr.r);

  print('A = Q * R');
  print(aqr.q.matmul(aqr.r));

  Double2D b = (a.cols[0].plus(a.cols[1] * 5)).toCol();
  print('B');
  print(b);

  print('x = A^-1 * b');
  print(aqr.solve(b));
}
