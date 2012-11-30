mata:

for (rvars = 1; rvars <= 3; rvars++) {
printf("\nRUNNING TEST CODE WITH " + strofreal(rvars) + " RANDOM VARIABLES\n");


/*
for (rv = 1; rv <= rvars; rv++) {
    printf("COL="+strofreal(rv)+"\n");
    for (rv2 = rv; rv2 <= rvars; rv2 = rv2 + 1 + max((rvars - rv, 0))) {
        printf("    ROW="+strofreal(rv2)+"\n");
        //printf("    increasing row by " + strofreal(rvars-rv) + "\n");
        //printf("replacing S2xx[|"+strofreal(rv2)+","+strofreal(rv)+"\"+strofreal(rv2+nRows-1)+","+strofreal(rv+nCols-1)+"|]\n");
    }
}
*/

printf("Test 1 -----------------\n");
iCol = 1
for (nCols = rvars; nCols >= 1; nCols--) {
    rv = rvars - nCols + 1
    iRow = iCol
    for (nRows = nCols; nRows >= 1; nRows--) {
        rv2 = rvars - nRows + 1
        printf("        replacing S2xx[|"+strofreal(iRow)+","+strofreal(iCol)+"\"+strofreal(iRow+nRows-1)+","+strofreal(iCol+nCols-1)+"|] with Rvar["+strofreal(rv2)+","+strofreal(rv)+"]\n");
        iRow = iRow + nRows
    }
    iCol = iCol + nCols
}

printf("Test 2 -----------------\n");
iCol = 1
/*
for (nCols = rvars; nCols >= 1; nCols--) {
    rv = rvars - nCols + 1
*/
for (rv = 1; rv <= rvars; rv++) {
    nCols = rvars - rv + 1
    iRow = iCol
    /*
    for (nRows = nCols; nRows >= 1; nRows--) {
        rv2 = rvars - nRows + 1
    */
    for (rv2 = rv; rv2 <= rvars; rv2++) {
        nRows = rvars - rv2 + 1
        printf("        replacing S2xx[|"+strofreal(iRow)+","+strofreal(iCol)+"\"+strofreal(iRow+nRows-1)+","+
                strofreal(iCol+nCols-1)+"|] with Rvar["+strofreal(rv2)+","+strofreal(rv)+"]\n");
        iRow = iRow + nRows
    }
    iCol = iCol + nCols
}

}

end

***