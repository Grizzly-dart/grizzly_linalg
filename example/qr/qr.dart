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

  final xQR = qr(x);
  print(xQR.q);
  print(xQR.r);

  print(xQR.q.shape);
  print(xQR.r.shape);
}