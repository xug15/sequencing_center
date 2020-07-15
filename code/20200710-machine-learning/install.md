# install.md
## Docker 安装
```sh
# Download tensorflow
docker pull tensorflow/tensorflow
docker pull tensorflow/tensorflow:latest-jupyter

# Start a CPU-only container
$ docker run -it --rm tensorflow/tensorflow bash

# Start a GPU container, using the Python interpreter.
$ docker run -it --rm --runtime=nvidia tensorflow/tensorflow:latest-gpu python

# Run a Jupyter notebook server with your own notebook directory (assumed here to be ~/notebooks). To use it, navigate to localhost:8888 in your browser.
$ docker run -it --rm -v $(realpath ~/notebooks):/tf/notebooks -p 6006:6006 -p 8888:8888 tensorflow/tensorflow:latest-jupyter
# xugang computer example.
docker run -it --rm -v /Users/xugang/Documents/c-pycharm/git/sequencing_center/code/20200710-machine-learning/notebooks:/tf/notebooks  -p 6006:6006 -p 8888:8888 tensorflow/tensorflow:latest-jupyter

```
在浏览器中打开：
```sh
http://localhost:8888/
http://localhost:6006/
```

## 练习
```sh
# cpu
docker pull cargo.caicloud.io/tensorflow/tensorflow:0.12.0
# run
docker run -it -v /Users/xugang/Documents/c-pycharm/git/sequencing_center/code/20200710-machine-learning/notebooks:/Deep_Learning_with_TensorFlow  -p 6006:6006 -p 8888:8888 cargo.caicloud.io/tensorflow/tensorflow:0.12.0
#gpu
docker pull cargo.caicloud.io/tensorflow/tensorflow:0.12.0-gpu
# run
nvidia-docker run -it -p 6006:6006 -p 8888:8888 cargo.caicloud.io/tensorflow/tensorflow:0.12.0-gpu
```

## 打开浏览器
密码可以在启动docker 的终端里面找到。
```sh
http://localhost:8888/
172.17.0.8:6006
```

