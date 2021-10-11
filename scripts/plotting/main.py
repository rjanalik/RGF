#!/usr/bin/env python3
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import argparse

sns.set()
# df = pd.read_csv("results.csv", sep="\t")
# df["speed_up"] = df["mean"][0]/df["mean"]
# df["std_speed_up"] = df["std"]/df["std"][0]

# plt.errorbar(df["num_threads"], df["mean"], yerr=df["std"])
# plt.savefig("runtimes.svg")

# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt.errorbar(df["num_threads"], df["speed_up"], yerr=df["std_speed_up"], label='Array reduction')
# x = np.linspace(0, df.shape[0], 1000)
# plt.plot(x, x, linestyle='dashed', label="Ideal-speed up")
# plt.title("16 bins [XeonGold_6150]")
# plt.ylabel("[threads / speed-up]", fontsize=14, rotation=0)
# ax.yaxis.set_label_coords(0.02, 1.10)
# plt.legend(loc='lower right')
# plt.savefig("speed_up.svg")


if __name__ == '__main__':

# Create results in here unless we specify a logdir
    # if FLAGS.logdir is not None and not os.path.exists(FLAGS.logdir):
    #     os.mkdir(FLAGS.logdir)

    x_var = 'no'
    df = pd.read_csv ('../../results/ghcn/results.csv', sep="\t")
    df = df.sort_values(["ns", "no"], ascending=False)
    df = df.loc[df["ns"]==4002]
    TYPE = "t"
    # df["speed_up"] = df["mean"][0]/df["mean"]
    # df["std_speed_up"] = df["std"]/df["std"][0]
    # df.loc[df["ns"]==40]
    fig = plt.figure()

    ax = fig.add_subplot(111)
    for version in df["RGF_Version"].unique():
        plt.scatter(df.loc[df['RGF_Version'] == version][x_var], df.loc[df['RGF_Version'] == version][TYPE+'_factorize'], label=version+' - factorize')
        # plt.scatter(df.loc[df['RGF_Version'] == version][x_var], df.loc[df['RGF_Version'] == version][TYPE+'_solve'], label=version+' - solve')
        # plt.scatter(df.loc[df['RGF_Version'] == version][x_var], df.loc[df['RGF_Version'] == version][TYPE+'_inv'], label=version+' - inv')
    # plt.scatter(df[x_var], df[TYPE+'_solve'], label='solve')
    # plt.scatter(df[x_var], df[TYPE+'_inv'], label='inversion')
    plt.ylabel("runtime [s]", fontsize=14, rotation=0)
    ax.yaxis.set_label_coords(0.02, 1.04)
    plt.xlabel("size [n]", fontsize=14, rotation=0)
    plt.legend(loc='lower right')
    plt.title("NVIDIA Quadro GV100")
    # x = np.linspace(0, df.shape[0], 1000)
    # plt.plot(x, x, linestyle='dashed', label="Ideal-speed up")
    plt.show()
