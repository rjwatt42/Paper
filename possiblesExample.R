########################

setWorld(TRUE,"Single","z",0.3)
setBrawEnv("RZ",'z')
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=NULL))
g<-showPossible(p,showType="Samples")
print(g)

#######################

setWorld(TRUE,"Single","z",0.3,pRplus=0.5)
setBrawEnv("RZ",'z')
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=NULL))
g<-showPossible(p,showType="Samples")
print(g)

#######################

setWorld(TRUE,"Exp","z",0.3,pRplus=1)
setBrawEnv("RZ",'z')
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=NULL))
g<-showPossible(p,showType="Samples")
print(g)

#######################

setWorld(TRUE,"Single","z",0.3,pRplus=0.5)
setBrawEnv("RZ",'z')
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=0.2))
g<-showPossible(p,showType="Samples")
print(g)

#######################

setWorld(TRUE,"Single","z",0.3,pRplus=0.25)
setBrawEnv("RZ",'z')
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=0.2))
g<-showPossible(p,showType="Samples")
print(g)

#######################
setWorld(TRUE,"Exp","z",0.2)
setDesign(sN=50,sNRand = FALSE)

p<-doPossible(makePossible(targetSample=0.5,sigOnly=FALSE))
g<-showPossible(p,showType="Samples",doTextResult = "PDF")
print(g)
p<-doPossible(makePossible(targetSample=0.5,sigOnly=TRUE))
g<-showPossible(p,showType="Samples",doTextResult = "PDF")
print(g)

