def getLifetimeX(mGluino, F, mX):
    return 1.e-6*pow(sqrt(F)/1000.e12,4)*pow(3.e12/mGluino,2)*pow(1.e9/mX,3)
