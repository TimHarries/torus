FROM gcc:7.2
WORKDIR /app
COPY . /app
RUN /app/buildtorus single cfitsio=no debug=yes
