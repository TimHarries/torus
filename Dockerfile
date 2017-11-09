FROM dockermpi
WORKDIR /app
COPY . /app
RUN /app/buildtorus mpi=yes cfitsio=no debug=yes 
