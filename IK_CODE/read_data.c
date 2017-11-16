/* This program is special procedure to read in the data of molecule as
in "data_mol". */

#include <stdio.h>

mygetc(fp1)
FILE* fp1;
{
char a;
a=getc(fp1);
while ((a != '(')&& ( a!=EOF))
{
a=getc(fp1);
}
return a;
}

mygetname(fp1)
FILE* fp1;
{
char g;
g=getc(fp1);
g=getc(fp1);
while ((g != 'c') && (g != 'n') && (g != EOF) )
{
g=getc(fp1);
}
putchar(g);
printf("\n");
return g;

}
