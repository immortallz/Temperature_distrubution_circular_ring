# Компилятор и флаги
CXX = g++
CXXFLAGS = -Wall -O2

# Определение команды для удаления файлов в зависимости от ОС
ifeq ($(OS),Windows_NT)
    RM = del /F /Q
else
    RM = rm -f
endif

# Исполняемый файл и объектные файлы
TARGET = a.exe
OBJS = elliptic.o main.o

# Цель по умолчанию
all: $(TARGET)

# Линковка
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(CXXFLAGS) -o $(TARGET)

# Компиляция файлов
elliptic.o: elliptic.cpp r.h
	$(CXX) $(CXXFLAGS) -c elliptic.cpp

main.o: main.cpp r.h
	$(CXX) $(CXXFLAGS) -c main.cpp

# Запуск программы с последующей очисткой временных файлов
run: $(TARGET)
	./$(TARGET)
	$(MAKE) clean

# Очистка временных файлов
clean:
	$(RM) $(OBJS) $(TARGET)

.PHONY: all run clean