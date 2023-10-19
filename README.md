
# Parallel PageRank algorithms  

## Table of Contents

- [Parallel PageRank algorithms](#parallel-pagerank-algorithms)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Project Overview](#project-overview)
  - [Building and Running the project](#building-and-running-the-project)


## Introduction

This project demonstrates parallel implementations of the PageRank algorithm using C++11 threads and provides a choice between Pull-Based and Push-Based strategies. Additionally, the project introduces four different task decomposition and mapping strategies.

## Project Overview

In this project, you will find three task decompositions and mapping strategies:  

1. **Vertex-based Decomposition and Static Mapping** 
2. **Edge-based Decomposition** 
3. **Vertex-based Decomposition & Dynamic Mapping**
4. **Vertex-based Decomposition & Coarse-Grained Dynamic Mapping**

## Building and Running the project

1. **Clone the repository** 
```bash 
$ git clone git@github.com:abay-kulamkadyr/parallel_PageRank_algorithms.git
```
2. **To Build**: Run `make` to build the project as follows: 
```bash 
$ make page_rank_<push OR pull>_parallel
```
3.  **Run**: 
```bash 
./page_rank_<pull OR Push>_parallel --nThreads <number_of_threads> --nIterations <number_of_iterations> --inputFile <path_to_input_graph> --strategy <strategy_number> --granularity <granularity_num>
```