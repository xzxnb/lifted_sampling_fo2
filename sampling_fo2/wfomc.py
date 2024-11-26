from __future__ import annotations
from collections import defaultdict

from enum import Enum
from functools import reduce
import math
import os
import argparse
import logging
import logzero

from logzero import logger
from typing import Callable
from contexttimer import Timer

from sampling_fo2.cell_graph.cell_graph import build_cell_graphs
from sampling_fo2.problems import WFOMCSProblem

from sampling_fo2.utils import MultinomialCoefficients, multinomial, \
    multinomial_less_than, RingElement, Rational, round_rational
from sampling_fo2.cell_graph import CellGraph, Cell
from sampling_fo2.context import WFOMCContext
from sampling_fo2.parser import parse_input
from sampling_fo2.fol.syntax import Const, Pred, QFFormula, PREDS_FOR_EXISTENTIAL
from sampling_fo2.utils.polynomial import coeff_dict, create_vars, expand

from parser.mln_parser import parse as mln_parse
from sampling_fo2.problems import WFOMCSProblem, MLN_to_WFOMC
import copy
class Algo(Enum):
    STANDARD = 'standard'
    FASTER = 'faster'
    FASTERv2 = 'fasterv2'
    INCREMENTAL = 'incremental'

    def __str__(self):
        return self.value

# 计算给定配置的权重，使用更快的算法
# config: 配置列表
# cell_weights: 单元权重列表
# edge_weights: 边权重矩阵
# 返回配置的权重
def get_config_weight_standard_faster(config: list[int],
                                      cell_weights: list[RingElement],
                                      edge_weights: list[list[RingElement]]) \
        -> RingElement:
    res = Rational(1, 1)
    for i, n_i in enumerate(config):
        if n_i == 0:
            continue
        n_i = Rational(n_i, 1)
        res *= cell_weights[i] ** n_i
        res *= edge_weights[i][i] ** (n_i * (n_i - 1) // Rational(2, 1))
        for j, n_j in enumerate(config):
            if j <= i:
                continue
            if n_j == 0:
                continue
            n_j = Rational(n_j, 1)
            res *= edge_weights[i][j] ** (n_i * n_j)
    return res

# 计算给定单元图和配置的权重，使用标准算法
# cell_graph: 单元图
# cell_config: 单元配置字典
# 返回配置的权重
def get_config_weight_standard(cell_graph: CellGraph,
                               cell_config: dict[Cell, int]) -> RingElement:
    res = Rational(1, 1)
    for i, (cell_i, n_i) in enumerate(cell_config.items()):
        if n_i == 0:
            continue
        res = res * cell_graph.get_cell_weight(cell_i) ** n_i
        res = res * cell_graph.get_two_table_weight(
            (cell_i, cell_i)
        ) ** (n_i * (n_i - 1) // 2)
        for j, (cell_j, n_j) in enumerate(cell_config.items()):
            if j <= i:
                continue
            if n_j == 0:
                continue
            res = res * cell_graph.get_two_table_weight(
                (cell_i, cell_j)
            ) ** (n_i * n_j)
    # logger.debug('Config weight: %s', res)
    return res

# 使用优化的单元图和更快的算法计算加权一阶模型计数（WFOMC）
# formula: 公式
# domain: 域
# get_weight: 获取权重的回调函数
# modified_cell_sysmmetry: 是否使用修改后的单元对称性
# 返回WFOMC结果
def faster_wfomc(formula: QFFormula,
                 domain: set[Const],
                 get_weight: Callable[[Pred], tuple[RingElement, RingElement]],
                 modified_cell_symmetry: bool = False) -> RingElement:
    domain_size = len(domain)
    res = Rational(0, 1)
    for opt_cell_graph, weight in build_cell_graphs(
        formula, get_weight, optimized=True,
        domain_size=domain_size,
        modified_cell_symmetry=modified_cell_symmetry
    ):
        cliques = opt_cell_graph.cliques
        nonind = opt_cell_graph.nonind
        i2_ind = opt_cell_graph.i2_ind
        nonind_map = opt_cell_graph.nonind_map

        res_ = Rational(0, 1)
        with Timer() as t:
            for partition in multinomial_less_than(len(nonind), domain_size):
                mu = tuple(partition)
                if sum(partition) < domain_size:
                    mu = mu + (domain_size - sum(partition),)
                coef = MultinomialCoefficients.coef(mu)
                body = Rational(1, 1)

                for i, clique1 in enumerate(cliques):
                    for j, clique2 in enumerate(cliques):
                        if i in nonind and j in nonind:
                            if i < j:
                                body = body * opt_cell_graph.get_two_table_weight(
                                    (clique1[0], clique2[0])
                                ) ** (partition[nonind_map[i]] *
                                      partition[nonind_map[j]])

                for l in nonind:
                    body = body * opt_cell_graph.get_J_term(
                        l, partition[nonind_map[l]]
                    )
                    if not modified_cell_symmetry:
                        body = body * opt_cell_graph.get_cell_weight(
                            cliques[l][0]
                        ) ** partition[nonind_map[l]]

                opt_cell_graph.setup_term_cache()
                mul = opt_cell_graph.get_term(len(i2_ind), 0, partition)
                res_ = res_ + coef * mul * body
        res = res + weight * res_
    logger.info('WFOMC time: %s', t.elapsed)
    return res

<<<<<<< HEAD
# 使用标准算法计算加权一阶模型计数（WFOMC）
# formula: 公式
# domain: 域
# get_weight: 获取权重的回调函数
# 返回WFOMC结果
=======

>>>>>>> 560c316bf926cb1018f83829f80c5cf7ba4c0328
def standard_wfomc(formula: QFFormula,
                   domain: set[Const],
                   get_weight: Callable[[Pred], tuple[RingElement, RingElement]]) -> RingElement:
    # cell_graph.show()
    res = Rational(0, 1)
    domain_size = len(domain)
    for cell_graph, weight in build_cell_graphs(formula, get_weight):
        res_ = Rational(0, 1)
        cells = cell_graph.get_cells()
        n_cells = len(cells)
        for partition in multinomial(n_cells, domain_size):
            coef = MultinomialCoefficients.coef(partition)
            cell_config = dict(zip(cells, partition))
            # logger.debug(
            #     '=' * 15 + ' Compute WFOMC for the partition %s ' + '=' * 15,
            #     dict(filter(lambda x: x[1] != 0, cell_config.items())
            # ))
            res_ = res_ + coef * get_config_weight_standard(
                cell_graph, cell_config
            )
        res = res + weight * res_
    return res


def incremental_wfomc(formula: QFFormula,
                      domain: set[Const],
                      get_weight: Callable[[Pred],
                                           tuple[RingElement, RingElement]],
                      leq_pred: Pred = None) -> RingElement:
    res = Rational(0, 1)
    domain_size = len(domain)
    for cell_graph, weight in build_cell_graphs(
            formula, get_weight, leq_pred=leq_pred
    ):
        # cell_graph.show()
        cells = cell_graph.get_cells()
        n_cells = len(cells)
        domain_size = len(domain)

        table = dict(
            (tuple(int(k == i) for k in range(n_cells)),
             cell_graph.get_cell_weight(cell))
            for i, cell in enumerate(cells)
        )
        for _ in range(domain_size - 1):
            old_table = table
            table = dict()
            for j, cell in enumerate(cells):
                w = cell_graph.get_cell_weight(cell)
                for ivec, w_old in old_table.items():
                    w_new = w_old * w * reduce(
                        lambda x, y: x * y,
                        (
                            cell_graph.get_two_table_weight((cell, cells[k]))
                            ** int(ivec[k]) for k in range(n_cells)
                        ),
                        Rational(1, 1)
                    )
                    ivec = list(ivec)
                    ivec[j] += 1
                    ivec = tuple(ivec)

                    w_new = w_new + table.get(ivec, Rational(0, 1))
                    table[tuple(ivec)] = w_new
        res = res + weight * sum(table.values())

    # if leq_pred is not None:
    #     res *= Rational(math.factorial(domain_size), 1)
    return res

# 预计算存在量词的权重
# cell_graph: 单元图
# domain_size: 域大小
# context: WFOMC上下文
# 返回存在量词权重的字典
def precompute_ext_weight(cell_graph: CellGraph, domain_size: int,
                          context: WFOMCContext) \
        -> dict[frozenset[tuple[Cell, frozenset[Pred], int]], RingElement]:
    existential_weights = defaultdict(lambda: Rational(0, 1))
    cells = cell_graph.get_cells()
    eu_configs = []
    for cell in cells:
        config = []
        for domain_pred, tseitin_preds in context.domain_to_evidence_preds.items():
            if cell.is_positive(domain_pred):
                config.append((
                    cell.drop_preds(
                        prefixes=PREDS_FOR_EXISTENTIAL), tseitin_preds
                ))
        eu_configs.append(config)

    cell_weights, edge_weights = cell_graph.get_all_weights()

    for partition in multinomial(len(cells), domain_size):
        # res = get_config_weight_standard(
        #     cell_graph, dict(zip(cells, partition))
        # )
        res = get_config_weight_standard_faster(
            partition, cell_weights, edge_weights
        )
        eu_config = defaultdict(lambda: 0)
        for idx, n in enumerate(partition):
            for config in eu_configs[idx]:
                eu_config[config] += n
        eu_config = dict(
            (k, v) for k, v in eu_config.items() if v > 0
        )
        existential_weights[
            frozenset((*k, v) for k, v in eu_config.items())
        ] += (Rational(MultinomialCoefficients.coef(partition), 1) * res)
    # remove duplications
    for eu_config in existential_weights.keys():
        dup_factor = Rational(MultinomialCoefficients.coef(
            tuple(c[2] for c in eu_config)
        ), 1)
        existential_weights[eu_config] /= dup_factor
    return existential_weights

<<<<<<< HEAD
# 根据指定算法计算加权一阶模型计数（WFOMC）
# context: WFOMC上下文
# algo: 算法类型（标准或更快）
# 返回WFOMC结果
def wfomc(context: WFOMCContext, algo: Algo = Algo.STANDARD) -> Rational:
=======

def wfomc(problem: WFOMCSProblem, algo: Algo = Algo.STANDARD) -> Rational:
    # both standard and faster WFOMCs need precomputation
    if algo == Algo.STANDARD or algo == Algo.FASTER or \
            algo == algo.FASTERv2:
        MultinomialCoefficients.setup(len(problem.domain))

    context = WFOMCContext(problem)
    leq_pred = Pred('LEQ', 2)
    if leq_pred in context.formula.preds():
        logger.info('Linear order axiom with the predicate LEQ is found')
        logger.info('Invoke incremental WFOMC')
        algo = Algo.INCREMENTAL
    else:
        leq_pred = None

    with Timer() as t:
        if algo == Algo.STANDARD:
            res = standard_wfomc(
                context.formula, context.domain, context.get_weight
            )
        elif algo == Algo.FASTER:
            res = faster_wfomc(
                context.formula, context.domain, context.get_weight
            )
        elif algo == Algo.FASTERv2:
            res = faster_wfomc(
                context.formula, context.domain, context.get_weight, True
            )
        elif algo == Algo.INCREMENTAL:
            res = incremental_wfomc(
                context.formula, context.domain,
                context.get_weight, leq_pred
            )
    res = context.decode_result(res)
    logger.info('WFOMC time: %s', t.elapsed)
    return res


def count_distribution(problem: WFOMCSProblem, preds: list[Pred],
                       algo: Algo = Algo.STANDARD) \
        -> dict[tuple[int, ...], Rational]:
    context = WFOMCContext(problem)
    # both standard and faster WFOMCs need precomputation
    if algo == Algo.STANDARD or algo == Algo.FASTER or \
            algo == algo.FASTERv2:
        MultinomialCoefficients.setup(len(problem.domain))
    leq_pred = Pred('LEQ', 2)
    if leq_pred in context.formula.preds():
        logger.info('Linear order axiom with the predicate LEQ is found')
        logger.info('Invoke incremental WFOMC')
        algo = Algo.INCREMENTAL
    else:
        leq_pred = None

    pred2weight = {}
    pred2sym = {}
    syms = create_vars('x0:{}'.format(len(preds)))
    for sym, pred in zip(syms, preds):
        if pred in pred2weight:
            continue
        weight = context.get_weight(pred)
        pred2weight[pred] = (weight[0] * sym, weight[1])
        pred2sym[pred] = sym
    context.weights.update(pred2weight)

>>>>>>> 560c316bf926cb1018f83829f80c5cf7ba4c0328
    if algo == Algo.STANDARD:
        res = standard_wfomc(
            context.formula, context.domain, context.get_weight
        )
    elif algo == Algo.FASTER:
        res = faster_wfomc(
            context.formula, context.domain, context.get_weight
        )
    elif algo == Algo.FASTERv2:
        res = faster_wfomc(
            context.formula, context.domain, context.get_weight, True
        )
<<<<<<< HEAD
    res = context.decode_result(res)
    return res

# 计算谓词的计数分布
# context: WFOMC上下文
# preds: 谓词列表
# algo: 算法类型（标准或更快）
# 返回计数分布的字典
def count_distribution(context: WFOMCContext, preds: list[Pred],
                       algo: Algo = Algo.STANDARD) \
        -> dict[tuple[int, ...], Rational]:
    pred2weight = {}
    pred2sym = {}
    syms = create_vars('x0:{}'.format(len(preds)))
    for sym, pred in zip(syms, preds):
        if pred in pred2weight:
            continue
        weight = context.get_weight(pred)
        pred2weight[pred] = (weight[0] * sym, weight[1])
        pred2sym[pred] = sym
    context.weights.update(pred2weight)
    if algo == Algo.STANDARD:
        res = standard_wfomc(
            context.formula, context.domain, context.get_weight
        )
    elif algo == Algo.FASTER:
        res = faster_wfomc(
            context.formula, context.domain, context.get_weight
=======
    elif algo == Algo.INCREMENTAL:
        res = incremental_wfomc(
            context.formula, context.domain,
            context.get_weight, leq_pred
>>>>>>> 560c316bf926cb1018f83829f80c5cf7ba4c0328
        )

    symbols = [pred2sym[pred] for pred in preds]
    count_dist = {}
    res = expand(res)
    for degrees, coef in coeff_dict(res, symbols):
        count_dist[degrees] = coef
    return count_dist


def parse_args():
    parser = argparse.ArgumentParser(
        description='WFOMC for MLN',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='mln file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default='./check-points')
    parser.add_argument('--algo', '-a', type=Algo,
                        choices=list(Algo), default=Algo.FASTER)
    parser.add_argument('--domain_recursive',
                        action='store_true', default=False,
                        help='use domain recursive algorithm '
                             '(only for existential quantified MLN)')
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()
    return args
def count_distribution_v2(context: WFOMCContext, preds1: list[Pred], preds2: list[Pred], mode: int,
                       algo: Algo = Algo.STANDARD) \
        -> dict[tuple[int, ...], Rational]:
    context_c = copy.deepcopy(context)
    pred2weight = {}
    pred2sym = {}
    preds = list(set(preds1+preds2))
    preds3 = {} #指定算哪部分的wmc
    if mode==1:
        preds3 = preds1
    else:
        preds3 = preds2
    syms = create_vars('x0:{}'.format(len(preds)))
    for sym, pred in zip(syms, preds):
        if pred in pred2weight:
            continue
        weight = context_c.get_weight(pred)
        if pred in preds3:
            pred2weight[pred] = (weight[0] * sym, 1)#False的weight赋1
        else:
            pred2weight[pred] = (sym, 1)      #preds2的先不管
        pred2sym[pred] = sym
    #把除preds即新preds以外的weight赋1
    for pred in context_c.formula.preds():
        if pred not in preds:
            pred2weight[pred] = (1,1)

    context_c.weights.update(pred2weight)
    aa = context_c.weights

    if algo == Algo.STANDARD:
        res = standard_wfomc(
            context_c.formula, context_c.domain, context_c.get_weight
        )
    elif algo == Algo.FASTER:
        res = faster_wfomc(
            context_c.formula, context_c.domain, context_c.get_weight
        )

    symbols = [pred2sym[pred] for pred in preds]
    count_dist = {}
    res = expand(res)
    for degrees, coef in coeff_dict(res, symbols):
        count_dist[degrees] = coef
    return count_dist
def MLN_TV(context: WFOMCContext, context1: WFOMCContext, preds1: list[Pred],
           context2: WFOMCContext, preds2: list[Pred]) -> Rational:
    #context1是mln1构造用来求Z1的formuler,preds1是新谓词,类似的context2，preds2，context是mln1、mln2合并构造的formuler

    Z1 = standard_wfomc(
        context1.formula, context1.domain, context1.get_weight
    )
    Z2 = standard_wfomc(
        context2.formula, context2.domain, context2.get_weight
    )
    count_dist1 = count_distribution_v2(context, list(preds1), list(preds2), 1)
    count_dist2 = count_distribution_v2(context, list(preds1), list(preds2), 2)
    res = Rational(0, 1)
    for key in count_dist1:
        res = res+abs(count_dist1[key]/Z1-count_dist2[key]/Z2)
    return res

def MLN_TV_V2(mln1: str,mln2: str) -> Rational:
    if mln1.endswith('.mln'):
        with open(mln1, 'r') as f:
            input_content = f.read()
        mln_problem1 = mln_parse(input_content)
        wfomcs_problem1 = MLN_to_WFOMC(mln_problem1,'@F')
        context1 = WFOMCContext(wfomcs_problem1)

    if mln2.endswith('.mln'):
        with open(mln2, 'r') as f:
            input_content = f.read()
        mln_problem2 = mln_parse(input_content)
        wfomcs_problem2 = MLN_to_WFOMC(mln_problem2, '@S')
        context2 = WFOMCContext(wfomcs_problem2)
    Z1 = standard_wfomc(
        context1.formula, context1.domain, context1.get_weight
    )
    Z2 = standard_wfomc(
        context2.formula, context2.domain, context2.get_weight
    )
    context = copy.deepcopy(context1)
    context.formula = context1.formula & context2.formula
    context.weights = {**context1.weights, **context2.weights}
    #这里没有对context的sentence做处理可能有隐患
    #利用context1.weights.keys只包含了新谓词，而context.formula.preds()包含所有谓词，可能有隐患
    #MLN_to_WFOMC貌似会对旧谓词的权重全赋为1
    a = context1.weights.keys()
    b = context2.weights.keys()
    count_dist1 = count_distribution_v2(context, list(a), list(b), 1)
    count_dist2 = count_distribution_v2(context, list(a), list(b), 2)
    res = Rational(0, 1)
    for key in count_dist1:
        res = res + abs(count_dist1[key] / Z1 - count_dist2[key] / Z2)
    return res






if __name__ == '__main__':
    # import sys
    # sys.setrecursionlimit(int(1e6))

<<<<<<< HEAD
    # args = parse_args()
    # if not os.path.exists(args.output_dir):
    #     os.makedirs(args.output_dir)
    # if args.debug:
    #     logzero.loglevel(logging.DEBUG)
    # else:
    #     logzero.loglevel(logging.INFO)
    # logzero.logfile('{}/log.txt'.format(args.output_dir), mode='w')
    #
    # with Timer() as t:
    #     problem = parse_input(args.input)
    # context = WFOMCContext(problem)
    # logger.info('Parse input: %ss', t)
    #
    # res = wfomc(
    #     context, algo=args.algo
    # )
    # logger.info('WFOMC (arbitrary precision): %s', res)
    # round_val = round_rational(res)
    # logger.info('WFOMC (round): %s (exp(%s))', round_val, round_val.ln())

    mln1 = "models\\friends-smokes.mln"
    mln2 = "models\\friends-smokes2.mln"
    res = MLN_TV_V2(mln1, mln2)
    print(res)



=======
    with Timer() as t:
        problem = parse_input(args.input)
    logger.info('Parse input: %ss', t)

    res = wfomc(
        problem, algo=args.algo
    )
    logger.info('WFOMC (arbitrary precision): %s', res)
    round_val = round_rational(res)
    logger.info('WFOMC (round): %s (exp(%s))', round_val, round_val.ln())
>>>>>>> 560c316bf926cb1018f83829f80c5cf7ba4c0328
