################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../dirent-1.20.1/examples/find.c \
../dirent-1.20.1/examples/locate.c \
../dirent-1.20.1/examples/ls.c \
../dirent-1.20.1/examples/updatedb.c 

OBJS += \
./dirent-1.20.1/examples/find.o \
./dirent-1.20.1/examples/locate.o \
./dirent-1.20.1/examples/ls.o \
./dirent-1.20.1/examples/updatedb.o 

C_DEPS += \
./dirent-1.20.1/examples/find.d \
./dirent-1.20.1/examples/locate.d \
./dirent-1.20.1/examples/ls.d \
./dirent-1.20.1/examples/updatedb.d 


# Each subdirectory must supply rules for building sources it contributes
dirent-1.20.1/examples/%.o: ../dirent-1.20.1/examples/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


