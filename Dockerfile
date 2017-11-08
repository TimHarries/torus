FROM gcc:7.2
WORKDIR /app
COPY . /app
ENV SYSTEM testsuite
RUN /app/buildtorus single cfitsio=no debug=yes
CMD ["/app/bin/torus.single"]
