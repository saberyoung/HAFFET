FROM ubuntu:20.04

ENV TZ=Asia/Shanghai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
WORKDIR /SNDA

ADD . /SNDA/

RUN sed -i 's/archive.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list \
    && sed -i 's/security.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list

RUN apt-get update -y \
    && apt-get install -y python3 python3-pip python3-dev \
    && apt-get install -y python3-tk

RUN cd megalut-sewpy \
    && python3 setup.py install \
    && cd ..

RUN pip3 install --upgrade pip -i https://pypi.tuna.tsinghua.edu.cn/simple \
    && pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple -r requirements.txt \
    && mkdir Data

RUN python3 setup.py install
