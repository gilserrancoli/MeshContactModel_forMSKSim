import casadi.*

x=MX.sym('x');
y1=(x^3)*5;
y2=0;
y=if_else(x>0,y1,y2)

fun=Function('f',{x},{y});

fun(0)

jac=jacobian(fun(x),x);
fjac=Function('fjac',{x},{jac});

C=CodeGenerator('gentest.c');
C.add(fun);
C.add(fun.jacobian);
C.generate();

