TYPE = 
FLAGS = -O3
OBJS = align84a.o utilpv34.o align44d.o align79b.o align73c.o align96e.o align71g.o align77f.o p_all28.o p_fun70a.o p_fun66b.o p_fun74c.o p_some12.o util25b.o util29a.o util11c.o heur14.o get_par2.o util18d.o util11e.o util8f.o util21g.o culeli3.o nodes11.o util18h.o util20i.o util4j.o util18k.o util19l.o util11m.o util13n.o util9o.o util4p.o fast8.o dalign24.o util2q.o utils.o bit-vect.o pvm.o util3r.o
LIBS = 
CC = gcc

malign$(TYPE) : $(OBJS)
	$(CC) $(FLAGS) $(OBJS) -o malign$(TYPE) $(LIBS)

align84a.o : align84a.c align3.h
	$(CC) $(FLAGS) -c $<

utilpv34.o : utilpv34.c align3.h
	$(CC) $(FLAGS) -c $<

align79b.o : align79b.c align3.h
	$(CC) $(FLAGS) -c $<

align73c.o : align73c.c align3.h
	$(CC) $(FLAGS) -c $<

align44d.o : align44d.c align3.h
	$(CC) $(FLAGS) -c $<

align96e.o : align96e.c align3.h
	$(CC) $(FLAGS) -c $<

align77f.o : align77f.c align3.h
	$(CC) $(FLAGS) -c $<

align71g.o : align71g.c align3.h
	$(CC) $(FLAGS) -c $<

p_all28.o : p_all28.c align3.h
	$(CC) $(FLAGS) -c $<

p_fun70a.o  : p_fun70a.c align3.h
	$(CC) $(FLAGS) -c $<

p_fun66b.o  : p_fun66b.c align3.h
	$(CC) $(FLAGS) -c $<

p_fun74c.o  : p_fun74c.c align3.h
	$(CC) $(FLAGS) -c $<

p_some12.o  : p_some12.c align3.h
	$(CC) $(FLAGS) -c $<

util25b.o  : util25b.c align3.h
	$(CC) $(FLAGS) -c $<

util29a.o  : util29a.c align3.h
	$(CC) $(FLAGS) -c $<

util11c.o   : util11c.c align3.h
	$(CC) $(FLAGS) -c $<

heur14.o    : heur14.c align3.h
	$(CC) $(FLAGS) -c $<

get_par2.o : get_par2.c align3.h
	$(CC) $(FLAGS) -c $<

util18d.o   : util18d.c align3.h
	$(CC) $(FLAGS) -c $<

util11e.o   : util11e.c align3.h
	$(CC) $(FLAGS) -c $<

util8f.o   : util8f.c align3.h
	$(CC) $(FLAGS) -c $<

util21g.o   : util21g.c align3.h
	$(CC) $(FLAGS) -c $<

culeli3.o   : culeli3.c align3.h
	$(CC) $(FLAGS) -c $<

nodes11.o   : nodes11.c align3.h
	$(CC) $(FLAGS) -c $<

util18h.o   : util18h.c align3.h
	$(CC) $(FLAGS) -c $<

util20i.o   : util20i.c align3.h
	$(CC) $(FLAGS) -c $<

util4j.o   : util4j.c align3.h
	$(CC) $(FLAGS) -c $<

util18k.o   : util18k.c align3.h
	$(CC) $(FLAGS) -c $<

util19l.o   : util19l.c align3.h
	$(CC) $(FLAGS) -c $<

util11m.o   : util11m.c align3.h
	$(CC) $(FLAGS) -c $<

util13n.o   : util13n.c align3.h
	$(CC) $(FLAGS) -c $<

util9o.o   : util9o.c align3.h
	$(CC) $(FLAGS) -c $<

util4p.o   : util4p.c align3.h
	$(CC) $(FLAGS) -c $<

fast8.o    : fast8.c align3.h utils.h yapp.h
	$(CC) $(FLAGS) -c $<

dalign24.o : dalign24.c align3.h
	$(CC) $(FLAGS) -c $<

util2q.o : util2q.c align3.h
	$(CC) $(FLAGS) -c $<

utils.o : utils.c align3.h utils.h
	$(CC) $(FLAGS) -c $<

bit-vect.o : bit-vect.c align3.h bit-vect.h
	$(CC) $(FLAGS) -c $<

pvm.o : pvm.c align3.h pvm.h
	$(CC) $(FLAGS) -c $<

util3r.o : util3r.c align3.h
	$(CC) $(FLAGS) -c $<

