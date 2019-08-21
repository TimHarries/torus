FROM tjharries/torusdocker:dockermpiv2
WORKDIR /app
COPY . /app
RUN /app/buildtorus mpi=yes cfitsio=no debug=yes SYSTEM=testsuite
