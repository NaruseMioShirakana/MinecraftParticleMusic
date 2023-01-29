# **我的世界粒子特效Midi播放器（以下简称MPMP）** <br><br>

## 介绍
---
### &emsp;项目名：我的世界粒子特效Midi播放器（粒子特效指令音乐）<br>
### &emsp;作者：训练有素的藤原白叶（藤原白葉）

### &emsp;测试视频地址：  [ **演示1** ](https://www.bilibili.com/video/BV11y4y1L7jT)      [ **演示2** ](https://www.bilibili.com/video/BV1uT4y1P7CX)

### &emsp;软件架构：
&emsp;&emsp;&emsp;&emsp;**1.**  仅包含.h头文件<br>
&emsp;&emsp;&emsp;&emsp;**2.**  需要自己创建C++（Cpp）文件并使用主函数调用头文件中函数<br>
&emsp;&emsp;&emsp;&emsp;**3.**  不存在.exe可执行文件（未来可能会加入）<br>

### &emsp;关于
&emsp;&emsp;&emsp;&emsp;MPMP是一个用于制作Minecraft原版指令特效粒子的工具。<br>
&emsp;&emsp;&emsp;&emsp;它可以通过MIDI的数据，自动确定坐标并生成一串mcfunction<br>
&emsp;&emsp;&emsp;&emsp;以在MC中显示粒子特效连线。<br>
&emsp;&emsp;&emsp;&emsp;MPMP可以简化制作MC特效音乐，省略制作过程中的重复步骤。<br> <br>

## 快速开始
---
### &emsp;安装：
#### &emsp;&emsp;1.  安装C++：
&emsp;&emsp;&emsp;&emsp; **i)：**   安装相关IDE，如&emsp;&emsp;&emsp;[DEVCPP(外部下载站，请注意计算机安全)](https://sourceforge.net/projects/orwelldevcpp/) &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;[VisualStudio(官网下载)](https://visualstudio.microsoft.com/zh-hans/)
#### &emsp;&emsp;2.  安装该项目：
&emsp;&emsp;&emsp;&emsp; **i）：** 点击本页面克隆/下载按钮，并点击下载ZIP<br>
&emsp;&emsp;&emsp;&emsp; **ii）：** 下载发行版中MC粒子特效Midi播放器并解压<br>
#### &emsp;&emsp;3.  配置MC环境：
&emsp;&emsp;&emsp;&emsp; **i) ：** 安装Minecraft - 1.16.5（Fabric）<br>
&emsp;&emsp;&emsp;&emsp; **ii) ：** 下载ColorBlock（1.16.5-Fabric）[点我下载](https://www.mcbbs.net/thread-917845-1-1.html)<br>
&emsp;&emsp;&emsp;&emsp; **iii) ：** 按照说明中的提示调用函数<br>

### &emsp;一些必须了解的基础知识：
#### &emsp;&emsp;1.  MC程序基础：
&emsp;&emsp;&emsp;&emsp; **i）：** `MC刻(tick)` MC游戏的时间单位为tick，它与现实中秒的对应关系为20:1，<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;即mc是以20的tps（tick per second）运行<br>
&emsp;&emsp;&emsp;&emsp; **ii）：** `MCFUNCTION` MC高版本提供的一个新功能，可以简单的理解为多个指令的<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;集合，它可以实现用一条指令执行多条指令<br>
&emsp;&emsp;&emsp;&emsp; **iii）：** `MCDATAPACK` MC高版本新功能，数据包，可以对原版内容进行修改，可以<br>&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;定义自己的Mcfunction<br>
#### &emsp;&emsp;2.  需要用到的MC指令：
&emsp;&emsp;&emsp;&emsp; **i）：**[execute](https://minecraft.fandom.com/zh/wiki/%E5%91%BD%E4%BB%A4/execute)<br>
&emsp;&emsp;&emsp;&emsp; **ii）：**[scorebroad](https://minecraft.fandom.com/zh/wiki/%E5%91%BD%E4%BB%A4/scoreboard)<br>
&emsp;&emsp;&emsp;&emsp; **iii）：**[particle](https://minecraft.fandom.com/zh/wiki/%E5%91%BD%E4%BB%A4/particle)<br>

### &emsp;使用说明：
  调用class，使用open函数打开写出文件的路径，然后调用class中的函数就可以了，参数分别是起点坐标和终点坐标。


## 将要做的
---
### &emsp;更多粒子特效：
#### &emsp;&emsp;基础线条：
&emsp;&emsp;&emsp;&emsp;<del>直线（普通直线和伪倒影效果直线）</del><br>
&emsp;&emsp;&emsp;&emsp;<del>抛物线（普通抛物线和伪倒影效果抛物线与结束点圈阵的抛物线）</del><br>
&emsp;&emsp;&emsp;&emsp;<del>正弦/余弦曲线。</del><br>
&emsp;&emsp;&emsp;&emsp;<del>螺线与螺线组。</del><br>
&emsp;&emsp;&emsp;&emsp;<del>平面连续相切圆和带有Y轴的连续相切圆</del><br>
&emsp;&emsp;&emsp;&emsp;<del>笛卡尔叶形线</del><br>
&emsp;&emsp;&emsp;&emsp;<del>星型线</del><br>
&emsp;&emsp;&emsp;&emsp;<del>倒映形线条</del><br>
&emsp;&emsp;&emsp;&emsp;<del>线条消失前散开</del><br>
&emsp;&emsp;&emsp;&emsp;贝塞尔曲线<br>
&emsp;&emsp;&emsp;&emsp;傅里叶变换<br>
&emsp;&emsp;&emsp;&emsp;钱学森弹道轨迹<br>
&emsp;&emsp;&emsp;&emsp;随机轨迹曲线<br>
&emsp;&emsp;&emsp;&emsp;<del>七影蝶飞行</del><br>
&emsp;&emsp;&emsp;&emsp;心形“抛物线”<br>
&emsp;&emsp;&emsp;&emsp;导弹作为线条头部<br>
#### &emsp;&emsp;音符特效：
&emsp;&emsp;&emsp;&emsp;<del>扩大的圆圈</del><br>
&emsp;&emsp;&emsp;&emsp;<del>螺旋线</del><br>
&emsp;&emsp;&emsp;&emsp;<del>星星</del><br>
&emsp;&emsp;&emsp;&emsp;<del>心形线</del><br>
&emsp;&emsp;&emsp;&emsp;<del>蝴蝶（翅膀会飞，会动的）</del><br>
&emsp;&emsp;&emsp;&emsp;米国和霓虹的城市名<br>
&emsp;&emsp;&emsp;&emsp;爆炸特效<br>
&emsp;&emsp;&emsp;&emsp;各种花里胡哨的特效<br>

### &emsp;MIDI：
&emsp;&emsp;&emsp;&emsp;MIDI文件读取<br>
&emsp;&emsp;&emsp;&emsp;MIDI文件编辑<br>
&emsp;&emsp;&emsp;&emsp;MIDI文件导出<br>
### &emsp;Playsound：
&emsp;&emsp;&emsp;&emsp;<del>自动处理播放声音指令</del><br>
&emsp;&emsp;&emsp;&emsp;适配各个网络上存在的音源资源包<br>
### &emsp;setblock：
&emsp;&emsp;&emsp;&emsp;<del>在每一个音符处放置方块</del><br>
### &emsp;可视化界面：
&emsp;&emsp;&emsp;&emsp;.exe可执行程序，可视化界面操作<br>
### &emsp;基岩版适配：
&emsp;&emsp;&emsp;&emsp;基岩版适配<br>
