set(SRCS
    daily.cc)

add_executable(daily
    ${SRCS})
target_link_libraries(daily
    sgp4)
