from pathlib import Path
import copy
from qlib.contrib.strategy.signal_strategy import WeightStrategyBase
from qlib.backtest.decision import TradeDecisionWO
from qlib.backtest.position import Position

full_ticker_list = [
    "ABBV",
    "ACN",
    "AEP",
    "AIZ",
    "ALLE",
    "AMAT",
    "AMP",
    "AMZN",
    "AVB",
    "AVY",
    "AXP",
    "BDX",
    "BF-B",
    "BMY",
    "BR",
    "CARR",
    "CDW",
    "CE",
    "CHTR",
    "CNC",
    "CNP",
    "COP",
    "CTAS",
    "CZR",
    "DG",
    "DPZ",
    "DRE",
    "DXC",
    "EWA",
    "EWC",
    "EWG",
    "EWH",
    "EWJ",
    "EWL",
    "EWQ",
    "EWT",
    "EWU",
    "EWY",
    "EWZ",
    "FTV",
    "GOOG",
    "GPC",
    "GSG",
    "HIG",
    "HIGH.L",
    "HST",
    "HYG",
    "IAU",
    "ICLN",
    "IEAA.L",
    "IEF",
    "IEFM.L",
    "IEMG",
    "IEUS",
    "IEVL.L",
    "IGF",
    "INDA",
    "IUMO.L",
    "IUVL.L",
    "IVV",
    "IWM",
    "IXN",
    "JPEA.L",
    "JPM",
    "KR",
    "LQD",
    "MCHI",
    "META",
    "MVEU.L",
    "OGN",
    "PG",
    "PPL",
    "PRU",
    "PYPL",
    "RE",
    "REET",
    "ROL",
    "ROST",
    "SEGA.L",
    "SHY",
    "SLV",
    "SPMV.L",
    "TLT",
    "UNH",
    "URI",
    "V",
    "VRSK",
    "VXX",
    "WRK",
    "XLB",
    "XLC",
    "XLE",
    "XLF",
    "XLI",
    "XLK",
    "XLP",
    "XLU",
    "XLV",
    "XLY",
    "XOM",
]


def path_gen(raw_path: str, dir: str = None, suffix: str = None) -> Path:
    path = Path(raw_path)

    if path.parent != Path("."):
        real_path = path.expanduser()
    else:
        if not dir:
            raise ValueError("Parameter dir is missing.")

        temp = Path.cwd()
        while True:
            if temp.match("M6"):
                break
            temp = temp.parent
        real_path = temp.joinpath(dir).joinpath(path)

    if real_path.suffix != suffix and suffix is not None:
        real_path = real_path.with_suffix(suffix)

    return real_path


class WeightProvidedStrategy(WeightStrategyBase):
    def _init_(
        self, *, hold_thresh=1, only_tradable=False, **kwargs,
    ):
        super().__init__(**kwargs)
        self.hold_thresh = hold_thresh
        self.only_tradable = only_tradable

    def generate_trade_decision(self, execute_result=None):
        # generate_trade_decision
        # generate_target_weight_position() and generate_order_list_from_target_weight_position() to generate order_list

        # get the number of trading step finished, trade_step can be [0, 1, 2, ..., trade_len - 1]
        trade_step = self.trade_calendar.get_trade_step()
        trade_start_time, trade_end_time = self.trade_calendar.get_step_time(trade_step)
        pred_start_time, pred_end_time = self.trade_calendar.get_step_time(trade_step, shift=1)
        pred_score = self.signal.get_signal(start_time=pred_start_time, end_time=pred_end_time)
        if pred_score is None:
            return TradeDecisionWO([], self)
        current_temp = copy.deepcopy(self.trade_position)
        assert isinstance(current_temp, Position)  # Avoid InfPosition
        target_weight = pred_score.to_dict()["score"]
        order_list = self.order_generator.generate_order_list_from_target_weight_position(
            current=current_temp,
            trade_exchange=self.trade_exchange,
            risk_degree=self.get_risk_degree(trade_step),
            target_weight_position=target_weight,
            pred_start_time=pred_start_time,
            pred_end_time=pred_end_time,
            trade_start_time=trade_start_time,
            trade_end_time=trade_end_time,
        )
        return TradeDecisionWO(order_list, self)
