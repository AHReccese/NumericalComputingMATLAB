Second Question
format long
coefficientA = input('enter the A coefficient: ');
coefficientB = input('enter the B coefficient: ');
coefficientC = input('enter the C coefficient: ');
bigNumber = 10 ^ 30;

if(coefficientA > bigNumber || coefficientB > bigNumber || coefficientC > bigNumber )
    coefficientA = coefficientA/bigNumber;
    coefficientB = coefficientB/bigNumber;
    coefficientC = coefficientC/bigNumber;
end
delta = coefficientB^2 - 4*coefficientA*coefficientC;
if(delta < 0 )
    disp("There are Complex roots!!!");
else
    
    if(coefficientA == 0 )
    firstRoot = ((-1)*coefficientC)/coefficientB;
    disp("It is not a quadratic equation it is linear Equation and the unique root is: ");
    disp(firstRoot);
    else    
    if(delta == 0)
    firstRoot = ((-1)*coefficientB)/(2*coefficientA);
    secondRoot = firstRoot;
    disp("There are repeated roots");
    disp("fisrtRoot = secondRoot: ");
    disp(firstRoot);
    end
    if(delta > 0 )
        deltaSqrt = sqrt(delta);
        if(coefficientB >=0)
            firstRoot = (2*coefficientC) / ((-1)*coefficientB -(deltaSqrt));
            secondRoot = ((-1)*coefficientB - deltaSqrt)/(2*coefficientA);
        else
            firstRoot = ((-1)*coefficientB + deltaSqrt)/(2*coefficientA);  
            secondRoot = (2*coefficientC)/((-1)*coefficientB +(deltaSqrt));
        end
        disp(" The firstRoot: " + firstRoot );
        disp(" The secondRoot: " + secondRoot);
    end
    end
    
end
firstExample: 
firstRoot = 0.5;     
   secondRoot = -1.3333;
secondExample:
firstRoot = 0.5;      
  secondRoot = -1.3333;
thirdExample: ( now it is a linear eqaution and there is just one root!)
unique root is = -1
forthExample:
firstRoot = 2.001;
secondRoot = 1.999;
fifthExample:
firstRoot = 9.999999999999999e+59 ;
secondRoot = 1;
Asked Question Answer
when the coefficientB is positive and more common when it is nearer to sqrt(delta) there is a possibleError of cancellation in (-b+sqrt(delta)) so in these types we use ( (2c)/(-b-sqrt(delta)) formula because as you see the chance of cancellation in the second formula becomes two weeker, and when the coefficientB is negative and again more common when its abstract value is nearer to sqrt(delta) there is a most probably cancellation error in the second formula ( 2c/(-b+sqrt(delta)))  so in this case we use the first formula (-b + sqrt(delta)) in order to avoid cancellation error.





