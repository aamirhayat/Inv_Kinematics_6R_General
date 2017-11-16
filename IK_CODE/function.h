#define coordinate double
#define DIM 25
#define PI 3.141592653
#define MAXG 24
#define MAXG2 48
#define MAXD 16
#define MAXD2 32

struct M_ENTRY {
    coordinate Cos;
    coordinate Sin;
    coordinate Const;
};

struct X_ENTRY {
    coordinate Quadratic;
    coordinate Linear;
    coordinate Const;
};
