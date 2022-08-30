FROM ubuntu:20.04 AS builder

RUN apt-get update && apt-get install -y make sudo

ADD . /build

WORKDIR /build
RUN make deps -j
RUN make grammars energies pairaligns -j

RUN virtualenv /knotify
RUN /knotify/bin/pip install -r /build/requirements.txt -r /build/wheel-requirements.txt
RUN /knotify/bin/pip install /build
RUN cp /build/*.so /knotify/lib/
RUN cp -r /build/pkenergy/hotknots/params /knotify/pkenergy-params

FROM ubuntu:20.04
LABEL maintainer="Angelos Kolaitis <neoaggelos@gmail.com>"

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 libmpfr6 libgomp1 libgsl23 \
    && rm -rf /var/lib/apt /var/cache/apt \
    && addgroup --gid 879 knotify \
    && adduser --system --uid 879 --group knotify

USER 879:879

ENV \
    KNOTIFY_YAEP_LIBRARY_PATH=/knotify/lib/libpseudoknot.so \
    KNOTIFY_BRUTEFORCE_LIBRARY_PATH=/knotify/lib/libbruteforce.so \
    KNOTIFY_SKIP_FINAL_AU_LIBRARY_PATH=/knotify/lib/libskipfinalau.so \
    KNOTIFY_CONSECUTIVE_PAIRALIGN_LIBRARY_PATH=/knotify/lib/libcpairalign.so \
    KNOTIFY_BULGES_LIBRARY_PATH=/knotify/lib/libbulges.so \
    KNOTIFY_PKENERGY=/knotify/lib/libpkenergy.so \
    KNOTIFY_PKENERGY_CONFIG_DIR=/knotify/pkenergy-params

COPY --from=builder /knotify /knotify
