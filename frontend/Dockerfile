FROM node:18-alpine
WORKDIR /app

# 依存関係
COPY package*.json ./
RUN npm install

# 全ファイルをコピー
COPY . .

# 念のため、実行時に public が存在するか確認するコマンドを挟む
CMD ["sh", "-c", "ls -R /app/public && npm start"]