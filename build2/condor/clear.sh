#!/bin/bash

if [ ! -d "./run" ]; then
  mkdir run
fi
rm ./run/*

if [ ! -d "./out" ]; then
  mkdir out
fi
rm ./out/*

if [ ! -d "./err" ]; then
  mkdir err
fi
rm ./err/*

if [ ! -d "./cmd" ]; then
  mkdir cmd
fi
rm ./cmd/*

if [ ! -d "./log" ]; then
  mkdir log
fi
rm ./log/*
