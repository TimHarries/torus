FROM tjharries/torusdocker:dockermpiv2
WORKDIR /app
COPY . /app
RUN /app/buildtorus mpi openmp hybrid single debug=yes SYSTEM=gfortran
