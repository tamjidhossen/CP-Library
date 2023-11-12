#!/bin/bash

if [ -z "$1" ]; then
  numLoop=1000
else
  numLoop=$1
fi

if [ -z "$2" ]; then
  doComp=1
else
  doComp=$2
fi

if [ "$doComp" -eq 1 ]; then
  echo "Compiling gen, myCode & correctCode..."

  g++ -std=c++20 gen.cpp -o gen
  g++ -std=c++20 myCode.cpp -o myCode
  g++ -std=c++20 correctCode.cpp -o correctCode

  echo "Done compiling."
fi

diff_found=""

for ((x=1; x<=$numLoop; x++)); do
  echo "Test Passed: $x"
  ./gen > generated_input
  ./myCode < generated_input > myAnswer
  ./correctCode < generated_input > correctAnswer

  # add -w option to ignore whitespace differences
  diff -w myAnswer correctAnswer > diagnostics
  if [ $? -ne 0 ]; then
    diff_found="y"
    break
  fi
done

if [ -n "$diff_found" ]; then
  echo ""
  echo "Difference Found"
  echo "Input:"
  cat generated_input
  echo -e "\nOutput:"
  cat myAnswer
  echo -e "\nExpected:"
  cat correctAnswer
else
  echo "All tests passed :D"
fi

rm generated_input
rm myAnswer
rm correctAnswer
