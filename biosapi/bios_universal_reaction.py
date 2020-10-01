

class BiosUniversalReaction:

    def __init__(self, json, api=None):
        self.json_data = json
        self.api = api
        self.database_reactions = None

    @property
    def reactions(self):
        if self.database_reactions is not None:
            return self.database_reactions

        ids = set()
        for m in self.json_data['members']:
            ids.add(m['id'])

        if self.api is not None:
            self.database_reactions = {}
            for r in self.api.get_reactions(ids):
                self.database_reactions[r.uid] = r
        return self.database_reactions