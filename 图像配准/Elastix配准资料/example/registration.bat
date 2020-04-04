@echo # elastix -f C:/Users/10446/Desktop/Registration2/elastix/example/exampleinput/fixed.mhd -m C:/Users/10446/Desktop/Registration2/elastix/example/exampleinput/moving.mhd -out C:/Users/10446/Desktop/Registration2/elastix/example/exampleoutput -p C:/Users/10446/Desktop/Registration2/elastix/example/exampleinput/parameters_Rigid.txt -p C:/Users/10446/Desktop/Registration2/elastix/example/exampleinput/parameters_BSpline.txt
@echo # C:/Users/10446/Desktop/Registration2/elastix/example/exampleinput
cd C:\Users\10446\Desktop\Convert3D\bin
c3d.exe -dicom-series-list-SeriesID "C:\Users\10446\Desktop\051_YE JIN SHENG_2661199_20180529_028Y_M\YE JIN SHENG_2661199_20180529_028Y_M_000\00000222.dcm"
pause
@echo # c3d.exe -dicom-series-read C:/Users/10446/Desktop/1  1.2.840.113619.2.55.3.2831164355.701.1526254298.770.3.31.250000512512 -type short -omc C:/Users/10446/Desktop/moving.mhd
@echo # c3d.exe -dicom-series-read C:/Users/10446/Desktop/2  1.2.392.200080.100.200.1597280977.14.11.15.211.253.0.10001.250000512512 -type short -omc C:/Users/10446/Desktop/fixed.mhd
 1.2.840.113619.2.55.3.2831164355.701.1526254298.770.3.31.250000512512