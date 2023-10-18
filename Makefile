ifdef USE_INT
MACRO = -DUSE_INT
endif

#compiler setup
CXX = g++
CXXFLAGS = -std=c++14 -O3 -pthread $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h core/graph.h core/quick_sort.h
SERIAL= page_rank_pull page_rank_push
PARALLEL= page_rank_pull_parallel page_rank_push_parallel page_rank_push_parallel_atomic
ALL= $(SERIAL) $(PARALLEL)


all : $(ALL)

% : %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
