################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../gcvspl.c 

CPP_SRCS += \
../Caller.cpp \
../DataSet.cpp \
../Globals.cpp \
../IntraSetRelation.cpp \
../Normalizer.cpp \
../Option.cpp \
../PercolatorCInterface.cpp \
../ResultHolder.cpp \
../SanityCheck.cpp \
../Scores.cpp \
../SetHandler.cpp \
../SqtSanityCheck.cpp \
../StdvNormalizer.cpp \
../UniNormalizer.cpp \
../main.cpp \
../ssl.cpp 

OBJS += \
./Caller.o \
./DataSet.o \
./Globals.o \
./IntraSetRelation.o \
./Normalizer.o \
./Option.o \
./PercolatorCInterface.o \
./ResultHolder.o \
./SanityCheck.o \
./Scores.o \
./SetHandler.o \
./SqtSanityCheck.o \
./StdvNormalizer.o \
./UniNormalizer.o \
./gcvspl.o \
./main.o \
./ssl.o 

C_DEPS += \
./gcvspl.d 

CPP_DEPS += \
./Caller.d \
./DataSet.d \
./Globals.d \
./IntraSetRelation.d \
./Normalizer.d \
./Option.d \
./PercolatorCInterface.d \
./ResultHolder.d \
./SanityCheck.d \
./Scores.d \
./SetHandler.d \
./SqtSanityCheck.d \
./StdvNormalizer.d \
./UniNormalizer.d \
./main.d \
./ssl.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


