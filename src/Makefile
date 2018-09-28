C    = cc  
OBJS =  menta.o fmax.o fmin.o allocate.o LocalDyna.o GlobalDyna.o  ReadInput.o OrderAlign.o \
        CreateMatrix4PAIR.o CreateMatrix4METHOD.o CreateMatrix4raw.o CreateMatrix4SF.o CreateMatrix4SQ.o \
        CreateMatrix4PCF.o CreateMatrix4PCQ.o CreateMatrix4LPCF.o CreateMatrix4LPCQ.o CreateMatrix4PROFSIM.o \
        CreateMatrix4PICASSO.o CreateMatrix4PICASSOQ.o CreateMatrix4COMPASS.o CreateMatrix4COMPASSW.o CreateMatrix4max.o\
        assign.o MaxData32.o MaxDataG32.o localtrace.o globaltrace.o WriteOutput.o extension.o aShell.o\
        iShell.o PrintIMatrix.o PrintFMatrix.o Seq2Profile.o PushSeq.o  PopSeq.o InitProfile.o ConsensusID.o\
        Align4Profile.o Align2Similarity.o Expand.o WriteProfile.o  MEntPSIC.o lambda.o NormalizePSSM.o PSSM2Matrix.o PrintPSSM.o\
        EVDmenta.o EVDisland.o statistic.o WriteProfileK.o SelectProfiles.o gshell.o
LIBS =  -lm 
.SUFFIXES: .c .o
.c.o:
	$(C) -c $*.c 
all: $(OBJS) 
	$(C) -o ../bin/menta $(OBJS) $(LIBS)
clean: 
	\rm *.o 
