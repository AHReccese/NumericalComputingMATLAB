Third Question

bisection
a = input('enter a: ')
b = input('enter b: ')
error = input('enter error: ')
while(true)
    c = (a+b)/2;
    if(abs((b-a)/2)<error || forthExample(c) == 0 )
        disp("root finded: ");
        disp(c);
        break;
    end
    if(forthExample(c)*forthExample(a) < 0)
        b = c;
    else
        a = c;
    end 
    disp("a: = ")
    disp(a);
    disp("b: =")
    disp(b);
    disp("c: =")
    disp(c);
end

Newton
x0 = input('enter the approximate root');
syms f(x)

% f(x) =  
% first = x^3-2*x-5;
% second = exp(-x)-x;
% third = xsin(x)-1;
% forth = x^3-3*x^2+3*x-1;

f(x) = x^3-3*x^2+3*x-1;
df = diff(f,x);
while(true)
    x1 = double(x0 - (f(x0)/double(df(x0))));
    disp("x(i): " + x1 );
    if(firstExample(x0) == 0 || abs(x1- x0)< 10^(-4))
        disp("root finded: ")
        disp(x1);
        break;
    end
    x0 = x1;
end

secant
x0 = input('enter first guess of root');
x1 = input('enter second guess of root');
while(true)
    rootFinder = x1 - (forthExample(x1)*(x1-x0))/(forthExample(x1)-forthExample(x0))
    if(abs(rootFinder - x1) < 10^(-4) )
        disp("root finded: ");
        disp(rootFinder);
        break
    else
        x0 = x1;
        x1 = rootFinder;
    end
end

Asked Question answers :
1.termination criterion: 
in the Bisection method as we read we know that when                                                              (abs(b-a)/2 < error ( given error))  -->  the root sequences are converged to the main root and we ignore the left of sequences.
in the Newton and Secant  Method i consider that when the abs(x(i+1)-x(i)) is less than a given error ( for example i consider these error is 10^-4)  the root sequences are converged to the result and we ignore the left sequences;

2.convergance rate of Newton & Secant & Bisection:
As we read in lesson we know that in Bisection Method  the root sequences converges to the result with the rate of convergance of  O((0.5)^n
again As we read in lesson we know that in Newton Method  the root sequences converges to the result with the rate of convergance of  1 + sqrt(2)
As we read in lesson we know that in Secant Method  the root sequences converges to the result with the rate of convergance of  (1 + sqrt(5))/2


Examples : 
function y = firstExample(x);
y = x^3 - 2*x - 5;
end
function y = secondExample(x);
y = exp(-x) - x;
end
function y = thirdExample(x);
y = x*sin(x)-1;
end
function y = forthExample(x);
y = x^3-3*x^2+3*x-1;
end
 
Bisection
Finding the Root of the Examples ( via Bisection ):
I consider that the acceptable error is given
first one :
a = 0
b = 3
error = 1.0000e-04
a: =    1.5000
b: =    3
c: =    1.5000
a: =    1.5000
b: =    2.2500
c: =    2.2500
a: =    1.8750
b: =    2.2500
c: =    1.8750
a: =    2.0625
b: =    2.2500
c: =    2.0625
a: =    2.0625
b: =    2.1563
c: =    2.1563
a: =    2.0625
b: =    2.1094
c: =    2.1094
a: =    2.0859
b: =    2.1094
c: =    2.0859
a: =    2.0859
b: =    2.0977
c: =    2.0977
a: =    2.0918
b: =    2.0977
c: =    2.0918
a: =    2.0918
b: =    2.0947
c: =    2.0947
a: =    2.0933
b: =    2.0947
c: =    2.0933
a: =    2.0940
b: =    2.0947
c: =    2.0940
a: =    2.0944
b: =    2.0947
c: =    2.0944
a: =    2.0945
b: =    2.0947
c: =    2.0945
root finded: 
    2.0946
Second One:
a = 0
b = 1
error = 1.0000e-04
a: =   0.5000
b: =     1
c: =    0.5000
a: =    0.5000
b: =    0.7500
c: =    0.7500
a: =     0.5000
b: =    0.6250
c: =   0.6250
a: =     0.5625
b: =    0.6250
c: =    0.5625
a: =     0.5625
b: =    0.5938
c: =    0.5938
a: =     0.5625
b: =    0.5781
c: =    0.5781
a: =     0.5625
b: =    0.5703
c: =    0.5703
a: =    0.5664
b: =   0.5703
c: =    0.5664
a: =     0.5664
b: =    0.5684
c: =    0.5684
a: =     0.5664
b: =    0.5674
c: =    0.5674
a: =     0.5669
b: =   0.5674
c: =   0.5669
a: =     0.5671
b: =   0.5674
c: =    0.5671
a: =     0.5671
b: =   0.5673
c: =    0.5673
root finded: 
    0.5672
Third One:
a = 0
b = 1.5708
error = 1.0000e-04
a: =     0.7854
b: =    1.5708
c: =    0.7854
a: =     0.7854
b: =    1.1781
c: =    1.1781
a: =     0.9817
b: =    1.1781
c: =    0.9817
a: =     1.0799
b: =    1.1781
c: =    1.0799
a: =     1.0799
b: =    1.1290
c: =    1.1290
a: =     1.1045
b: =    1.1290
c: =    1.1045
a: =     1.1045
b: =    1.1167
c: =    1.1167
a: =     1.1106
b: =   1.1167
c: =    1.1106
a: =     1.1137
b: =    1.1167
c: =    1.1137
a: =     1.1137
b: =   1.1152
c: =    1.1152
a: =     1.1137
b: =    1.1144
c: =    1.1144
a: =    1.1141
b: =    1.1144
c: =   1.1141
a: =     1.1141
b: =   1.1142
c: =   1.1142
root finded: 
    1.1141
Forth One:
a = 0
b = 3
error = 1.0000e-04
a: =     0
b: =   1.5000
c: =    1.5000
a: =    0.7500
b: =    1.5000
c: =    0.7500
a: =     0.7500
b: =    1.1250
c: =   1.1250
a: =     0.9375
b: =   1.1250
c: =   0.9375
a: =    0.9375
b: =    1.0313
c: =    1.0313
a: =    0.9844
b: =    1.0313
c: =    0.9844
a: =     0.9844
b: =    1.0078
c: =    1.0078
a: =    0.9961
b: =    1.0078
c: =    0.9961
a: =    0.9961
b: =    1.0020
c: =    1.0020
a: =     0.9990
b: =    1.0020
c: =    0.9990
a: =     0.9990
b: =    1.0005
c: =    1.0005
a: =     0.9998
b: =    1.0005
c: =   0.9998
a: =     0.9998
b: =    1.0001
c: =    1.0001
a: =     0.9999
b: =    1.0001
c: =    0.9999
root finded: 
    1.0000
secant
Finding the Root of the Examples ( via secant ):
I consider the maximum acceptable error is :x(i+1)-x(i) < 10^(-4)
first one : 
rootFinder = 2.5000
rootFinder = 2.0755
rootFinder = 2.0908
rootFinder = 2.0946
rootFinder = 2.0946
root finded: 
    2.0946
second one :
rootFinder = 0.6127
rootFinder = 0.5638
rootFinder = 0.5672
rootFinder = 0.5671
root finded: 
    0.5671
Third one :
rootFinder = 1
rootFinder = 1.1241
rootFinder = 1.1142
rootFinder = 1.1142
root finded: 
    1.1142
Forth one :
rootFinder = 0.3333
rootFinder = 0.4286
rootFinder = 0.5906
rootFinder = 0.6848
rootFinder = 0.7639
rootFinder = 0.8212
rootFinder = 0.8652
rootFinder = 0.8982
rootFinder = 0.9232
rootFinder = 0.9420
rootFinder = 0.9562
rootFinder = 0.9669
rootFinder = 0.9751
rootFinder = 0.9812
rootFinder = 0.9858
rootFinder = 0.9893
rootFinder = 0.9919
rootFinder = 0.9939
rootFinder = 0.9954
rootFinder = 0.9965
rootFinder = 0.9974
rootFinder = 0.9980
rootFinder = 0.9985
rootFinder = 0.9989
rootFinder = 0.9991
rootFinder = 0.9994
rootFinder = 0.9995
rootFinder = 0.9996
rootFinder = 0.9997
root finded: 
    0.9997
Newton
Finding the Root of the Examples ( via Newton ):
I consider the maximum acceptable error is :x(i+1)-x(i) < 10^(-4)
firstExample : x0 = 1.5;
x(i): 2.4737
x(i): 2.1564
x(i): 2.0966
x(i): 2.0946
x(i): 2.0946
root finded: 
    2.0946
secondExample : x0 = 0.5
x(i): 0.56631
x(i): 0.56714
x(i): 0.56714
root finded: 
    0.5671
ThirdExample : x0 = 1.4 
x(i): 1.0897
x(i): 1.1141
x(i): 1.1142
root finded: 
    1.1142
ForthExample : x0 = 2;
x(i): 1.6667
x(i): 1.4444
x(i): 1.2963
x(i): 1.1975
x(i): 1.1317
x(i): 1.0878
x(i): 1.0585
x(i): 1.039
x(i): 1.026
x(i): 1.0173
x(i): 1.0116
x(i): 1.0077
x(i): 1.0051
x(i): 1.0034
x(i): 1.0023
x(i): 1.0015
x(i): 1.001
x(i): 1.0007
x(i): 1.0005
x(i): 1.0003
x(i): 1.0002
x(i): 1.0001
root finded: 
    1.0001








