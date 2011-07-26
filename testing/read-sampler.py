import random
import sys

random.seed(0)

percent = float(sys.argv[1])
file = sys.argv[2]

fastqfile = open(file, 'r')
for line in fastqfile:
    if random.random() <= percent:
        sys.stdout.write(line)
        sys.stdout.write(fastqfile.next())
        sys.stdout.write(fastqfile.next())
        sys.stdout.write(fastqfile.next())
    else:
        fastqfile.next()
        fastqfile.next()
        fastqfile.next()


