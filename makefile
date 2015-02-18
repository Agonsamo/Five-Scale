# 5scale.make

# make file to compile 5scale program

# first name a variable objects for all object files
objects = 5-scale.o 4-scale3.o  distance.o liberty.o neyman.o p_gap_ax_b.o \
  ptg_sub.o  spheroid_ti.o equation1.o   ls.o      p_gap_branch.o \
  q1_ms.o    tree_size.o cone_gs.o   fo.o multiple_scattering.o   \
  poisson.o  q.o  triangle.o cone_ta.o   Initialise.o  overlap.o   pq_sub.o\
  spheroid_svg.o  xi.o cone_ti.o    nadir.o  parameters.o  ps.o\
  spheroid_ta.o

# name a variable sources for all source files
sources = 5-scale.c 4-scale3.c  distance.c liberty.c neyman.c p_gap_ax_b.c \
  ptg_sub.c  spheroid_ti.c equation1.c   ls.c      p_gap_branch.c \
  q1_ms.c    tree_size.c cone_gs.c   fo.c multiple_scattering.c   \
  poisson.c  q.c  triangle.c cone_ta.c   Initialise.c  overlap.c   pq_sub.c\
  spheroid_svg.c  xi.c cone_ti.c      nadir.c  parameters.c  ps.c\
  spheroid_ta.c


# now give target as makebeps with objects as variable dependencies + command line

5S: $(objects) 
	gcc -o 5S -g -lm $(objects)
	#g++ -o 5S -lm $(objects)

# list the dependencies for object files - those header files which help build objects

#5-scale.o					:direct.h
5-scale.o					:data.h
4-scale3.o				:data.h
tree_size.o				:data.h
spheroid_ti.o			:data.h
spheroid_ta.o			:data.h
spheroid_svg.o		:data.h
q1_ms.o						:data.h
q.o								:data.h
ptg_sub.o					:data.h
ps.o							:data.h
pq_sub.o					:data.h
poisson.o					:data.h
parameters.o			:data.h
p_gap_branch.o		:data.h
p_gap_ax_b.o			:data.h
overlap.o					:data.h
optical.o					:data.h	
neyman.o					:data.h
nadir.o						:data.h
multiple_scattering.o :data.h
ls.o							:data.h
liberty.o					:data.h
liberty.o					:liberty.h
input.o						:data.h
initialise.o			:data.h
fo.o							:data.h
distance.o				:data.h
cone_ti.o					:data.h
cone_ta.o					:data.h
cone_gs.o					:data.h

# how to build all object files from all dependent source files

$(objects): $(sources)
	gcc -c -g $(sources)
#gcc -c -g $(sources)
#	g++ -c -g $(sources)
# -g for debug -O for optimization

clean:
	rm $(objects) 5S

