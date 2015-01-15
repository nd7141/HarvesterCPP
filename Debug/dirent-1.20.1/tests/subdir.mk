################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../dirent-1.20.1/tests/t-dirent.c 

OBJS += \
./dirent-1.20.1/tests/t-dirent.o 

C_DEPS += \
./dirent-1.20.1/tests/t-dirent.d 


# Each subdirectory must supply rules for building sources it contributes
dirent-1.20.1/tests/%.o: ../dirent-1.20.1/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


