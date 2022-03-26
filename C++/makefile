EXP = expansion.o gaulag.o odeint.o spline.o
BBN = $(EXP) BBN.o nuclear.o stifbs.o ludcmp.o

fig1: fig1.o $(EXP)
	g++  fig1.o $(EXP)
fig2: fig2.o $(EXP)
	g++  fig2.o $(EXP)
fig4: fig4.o nuclear.o
	g++ fig4.o nuclear.o
fig5-6: fig5-6.o $(BBN)
	g++ fig5-6.o $(BBN)
fig7: fig7.o $(BBN)
	g++ fig7.o $(BBN)
fig8: fig8.o $(BBN)
	g++ fig8.o $(BBN)
