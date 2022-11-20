set(SRCS
    passpredict.cc)

add_executable(passpredict-sun
    ${SRCS})
target_link_libraries(passpredict-sun
    sgp4)
