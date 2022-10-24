import 'package:grizzly/grizzly.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  final a = [
    [1.0, 4.0],
    [2.0, 5.0],
  ];

  final b = a.matmul([5, 2].toCol());
  print(b);

  final LU xlu = lu(a);

  print(xlu.p);

  print(xlu.l);

  print(xlu.u);

  print(xlu.p.matmul(xlu.l).matmul(xlu.u));

  print(xlu.solve(b));
}
