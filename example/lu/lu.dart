import 'package:grizzly_array/grizzly_array.dart';
import 'package:grizzly_linalg/grizzly_linalg.dart';

main() {
  final a = new Double2D([
    [1.0, 4.0],
    [2.0, 5.0],
  ]);

  final b = a.matmul(new Double2D.aCol(<int>[5, 2]));
  print(b);

  final LU xlu = lu(a);

  print(xlu.p);

  print(xlu.l);

  print(xlu.u);

  print(xlu.p.matmul(xlu.l).matmul(xlu.u));

  print(xlu.solve(b));
}