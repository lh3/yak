CFLAGS=		-g -Wall -O2
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o bbf.o htab.o bseq.o misc.o sys.o 6gjdn.o \
			count.o qv.o triobin.o inspect.o
PROG=		yak
LIBS=		-lm -lz -lpthread

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

yak:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bbf.o: yak.h
bseq.o: bseq.h kseq.h
count.o: kthread.h yak-priv.h yak.h kseq.h
htab.o: kthread.h yak-priv.h yak.h khashl.h
inspect.o: ketopt.h yak-priv.h yak.h
kthread.o: kthread.h
main.o: ketopt.h yak-priv.h yak.h
misc.o: yak-priv.h yak.h
qv.o: kthread.h yak-priv.h yak.h bseq.h
sys.o: yak-priv.h yak.h
triobin.o: kthread.h ketopt.h bseq.h yak-priv.h yak.h
